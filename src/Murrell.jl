
module Murrell
"""
    A set of functions used to process data from the 'Murrell'
    gas medium triaxial, RIPL, UCL
"""

export M_file, M_read, H_mod, SFP_read, JR
using TDMSReader, XLSX, NoiseRobustDifferentiation
using DelimitedFiles

"""
    M_file(P, col::int64)

    Generate filepath and extract indices from spreadsheet log.xlsx

    #Arguments
    - 'P' : dictionary information pertaining to the 'Murrell'
    - 'col::Int64' : column of file of interest in accompanying spread sheet
"""
function M_file(P,col::Int64)
    if Sys.iswindows()
        xf = XLSX.readxlsx("C:\\Users\\cwaha\\Dropbox\\My PC (DESKTOP-8JF2H49)\\Documents\\UCL\\raw_lab data\\log.xlsx")
    else
        xf = XLSX.readxlsx("/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/raw_lab data/log.xlsx")
    end
    sh = xf[XLSX.sheetnames(xf)[1]]
    info = sh[:]
    P[:I]=[info[col,4],info[col,5],info[col,6]] # Indice at HP
    P[:P_exp] = info[col,2]
    P[:T] = info[col,7] # Temperature of experiment
    P[:εr] = info[col,8] # Strain rate of test
    P[:K_mm_kN] = info[col,9] # Machine compliance [mm kN^-1]
    P[:L_mm] = info[col,10] # Sample length
    P[:d_mm] = info[col,11] # Diameter of sample [m]
    fil = info[col,1]
    if Sys.iswindows()
        P[:fid] = "C:\\Users\\cwaha\\Dropbox\\My PC (DESKTOP-8JF2H49)\\Documents\\UCL\\raw_lab data\\"*fil*".tdms"
    else
        P[:fid] = "/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/raw_lab data/"*fil*".tdms"
    end
end
"""
    M_read(P)

    Open and read contents of .tdms file n.b. must be called after M_file
    #Arguments
    - 'P' : dictionary information pertaining to the 'Murrell'
"""
function M_read(P)
    #Read data from .tdms file
    tdmsIN = readtdms(P[:fid])
    P[:t_s] = tdmsIN.groups["Numeric"]["TimeStamp"].data
    P[:F_kN] = tdmsIN.groups["Numeric"]["Load"].data
    P[:Ua_mm] = tdmsIN.groups["Numeric"]["Displacement"].data
    P[:U1_mm] = tdmsIN.groups["Numeric"]["LVDT 1"].data
    P[:U2_mm] = tdmsIN.groups["Numeric"]["LVDT 2"].data
    P[:Pc1_MPa] = tdmsIN.groups["Numeric"]["Pc 700MPa"].data
    P[:Pc2_MPa] = tdmsIN.groups["Numeric"]["PC 1400 MPa"].data
    P[:Pf_MPa] = tdmsIN.groups["Numeric"]["Pore Pressure"].data
    P[:PpVol_mm3] = tdmsIN.groups["Numeric"]["PP vol"].data

    # Perform some data reduction
    P[:t_s_c] = P[:t_s].-P[:t_s][P[:I][1]]
    P[:F_kN_c] = P[:F_kN] .-P[:F_kN][P[:I][1]]
    P[:U_mm_c] = ((P[:U1_mm].+P[:U2_mm]).-(P[:U1_mm][P[:I][1]]+P[:U2_mm][P[:I][1]]))./2
    P[:U_mm_fc] = P[:U_mm_c] .-(P[:F_kN_c]*P[:K_mm_kN])
    P[:ε] = P[:U_mm_fc]./P[:L_mm]
    P[:Jr] = JR(P)
    P[:F_kN_j] = P[:F_kN_c] .-JR(P)
    P[:σ_MPa] = P[:F_kN_c]./(0.25e-6π*P[:d_mm]^2) .*1e-3
    P[:σ_MPa_j] = P[:F_kN_j]./(0.25e-6π*P[:d_mm]^2) .*1e-3
end
"""
    H_mod(P,ds::Int64)

    Calculate tangent moduli and yield stress using robust numerical differentiation N.B. must be called after M_read

    #Arguments
    - 'P' : dictionary information pertaining to the 'Murrell'
"""
function H_mod(P, ds::Int64)
    ε1 = Array{Float64,1}(undef,1)
    σ1 = Array{Float64,1}(undef,1)
    t1 = Array{Float64,1}(undef,1)
    ε = P[:ε][P[:I][1]:P[:I][2]]
    σ = P[:σ_MPa][P[:I][1]:P[:I][2]]
    t = P[:t_s][P[:I][1]:P[:I][2]]
    for i in 1:Int(floor(length(P[:ε][P[:I][1]:P[:I][2]])/ds))-1
        push!(ε1,sum(ε[Int(i*ds-ds+1):Int(i*ds)])/ds)
        push!(σ1,sum(σ[Int(i*ds-ds+1):Int(i*ds)])/ds)
        push!(t1,t[Int(i*ds)+Int(ceil(ds/2))])
    end
    deleteat!(ε1,1)
    deleteat!(σ1,1)
    deleteat!(t1,1)
    # dσ = differentiate(t1,σ1,TotalVariation(),0.4,5e-3,maxit=2000)
    dσ = TVRegDiff(σ1,100, 0.2, dx=P[:t_s][P[:I][1]+1]-P[:t_s][P[:I][1]],ε=1e-9)
    # dε = differentiate(t1,ε1,TotalVariation(),0.4,5e-3,maxit=2000)
    dε = TVRegDiff(ε1,100, 0.2, dx=P[:t_s][P[:I][1]+1]-P[:t_s][P[:I][1]],ε=1e-9)
    III = findfirst(σ1 .> 10) #exclude portion of experiment <0.2% strain, can give misleading results
    P[:h] = dσ[III:end]./dε[III:end]*1e-3 # Calculate the tangent modulus in GPa, excluding initial non-linearity
    P[:ε_h] = ε1
    h_max = maximum(P[:h][III:end])
    II = findfirst(P[:h] .== h_max)
    P[:E_GPa] = sum(P[:h][P[:h].>0.95h_max])/sum(P[:h].>0.95h_max) #youngs modulus about maximum hardening modulus
    P[:H_GPa] = sum(P[:h][P[:h].<0.4P[:E_GPa]])/sum(P[:h].<0.4P[:E_GPa])
    P[:ε_E] = ε1[II]
    P[:σ_E] = σ1[II]
    σ2 = σ1[III:end]
    ε2 = ε1[III:end]
    h = P[:h][III:end]
    P[:YSa_MPa] = σ2[findfirst(h.<0.8P[:E_GPa])] #80% yield
    P[:εa] = ε2[findfirst(h.<0.8P[:E_GPa])]
    P[:YSb_MPa] = σ2[findfirst(h.<0.5P[:E_GPa])] #50% yield
    P[:εb] = ε2[findfirst(h.<0.5P[:E_GPa])]
end
"""
    SFP_read()

    Read and return data from program SFP_furnace_control_V2.vi

    #Arguments
    - 'SFP' : empty dictionary to store extracted data
    - 'fil::string' : path to file
"""
function SFP_read(SFP,fil)
    path="/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/Furnace_calibration/SFP_logging/"*fil
    dat = readdlm(path)
    headers = dat[:,1]
    I = findall(x->x==headers[1],dat[:,1])
    SFP[:t] = dat[:,1]
    SFP[:PF] = dat[:,3]
    SFP[:S1] = dat[:,4]
    SFP[:S2] = dat[:,5]
    SFP[:S3] = dat[:,6]
    SFP[:TC1] = dat[:,7]
    SFP[:TC2] = dat[:,8]
    close(dat)
    deleteat!(SFP[:t],I)
    deleteat!(SFP[:PF],I)
    deleteat!(SFP[:S1],I)
    deleteat!(SFP[:S2],I)
    deleteat!(SFP[:S3],I)
    deletat!(SFP[:TC1],I)
    deletat!(SFP[:TC2],I)
end
"""
    JR(P)

    Calculate strength contribution of copper jacket for corrections

    #Arguments
    - 'P' : dictionary information pertaining to the 'Murrell'
"""
function JR(P)
    #Inputs: passed as a dictionary
    #T.......Temperature (°C)
    #ε.......Sample strain (-)
    #ε̇.......Strain rate during test (s⁻¹)
    #d.......Sample diameter (m)
    #L.......Sample length (m)
    EaJ = 197000 # Copper activation enthalpy
    n1J = 4.8 # Empirical factor 1
    n2J = 0.22 # Empirical factor 2
    ε̇j = P[:εr]*(P[:d_mm]/P[:L_mm]) # Jacket strain rate
    Jr = 2*10^((log10(ε̇j*exp(EaJ/(8.3145*(P[:T]+278)))))/n1J-n2J) # Copper flow stress at experiment conditions
    Ja = π*P[:d_mm]*P[:L_mm]*1e-6 # Jacket area
    JR = Jr.*Ja.*((1 .+P[:ε])*(P[:d_mm]/P[:L_mm])).*1e-3 # Force due to jacket assuming linear increase due to incremental strain
end
end
