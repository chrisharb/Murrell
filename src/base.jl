"""
    M_file(P, col::int64)

    Generate filepath and extract indices from spreadsheet log.xlsx

    Arguments
    - 'col::Int64' : column of file of interest in accompanying spread sheet
"""
function M_file(col::Int64)
    if Sys.iswindows()
        xf = XLSX.readxlsx("C:\\Users\\cwaha\\Dropbox\\My PC (DESKTOP-8JF2H49)\\Documents\\UCL\\raw_lab data\\log.xlsx")
    else
        xf = XLSX.readxlsx("/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/raw_lab data/log.xlsx")
    end
    sh = xf[XLSX.sheetnames(xf)[1]]
    info = sh[:]
    exp_info = Dict()
    exp_info[:I]=[info[col,4],info[col,5],info[col,6]] # Indice at HP
    exp_info[:P_exp] = info[col,2]
    exp_info[:T] = info[col,7] # Temperature of experiment
    exp_info[:εr] = info[col,8] # Strain rate of test
    exp_info[:K_mm_kN] = info[col,9] # Machine compliance [mm kN^-1]
    exp_info[:L_mm] = info[col,10] # Sample length
    exp_info[:d_mm] = info[col,11] # Diameter of sample [m]
    fil = info[col,1]
    if Sys.iswindows()
        fid = "C:\\Users\\cwaha\\Dropbox\\My PC (DESKTOP-8JF2H49)\\Documents\\UCL\\raw_lab data\\"*fil*".tdms"
    else
        fid = "/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/raw_lab data/"*fil*".tdms"
    end
    return fid, exp_info
end
"""
    M_read(P)

    Open and read contents of .tdms file
    #Arguments
    - fid : file path string
    return P
"""
function M_read(fid::String)
    tdmsIN = readtdms(fid)
    P = Dict()
    P[:t_s] = tdmsIN.groups["Numeric"]["TimeStamp"].data
    P[:F_kN] = tdmsIN.groups["Numeric"]["Load"].data
    P[:Ua_mm] = tdmsIN.groups["Numeric"]["Displacement"].data
    P[:U1_mm] = tdmsIN.groups["Numeric"]["LVDT 1"].data
    P[:U2_mm] = tdmsIN.groups["Numeric"]["LVDT 2"].data
    P[:Pc1_MPa] = tdmsIN.groups["Numeric"]["Pc 700MPa"].data
    P[:Pc2_MPa] = tdmsIN.groups["Numeric"]["PC 1400 MPa"].data
    P[:Pf_MPa] = tdmsIN.groups["Numeric"]["Pore Pressure"].data
    P[:PpVol_mm3] = tdmsIN.groups["Numeric"]["PP vol"].data
    return P
end

function t_conv(P)
    P[:t_date_time] = unix2datetime(datetime2unix(DateTime(1904,1,1,0,0,0)) .+Int.(Round.(P[:t_s]*1000)))
end

"""
    M_reduce!(P, exp_info)

    Reduce Murrell data
    #Arguments
    * P : dictionary containing raw data from the Murrell, must include indices of hit point and sample information
    * exp_info : dictionary containing experimental parameters
"""

function M_reduce!(P,exp_info)
    I1 = exp_info[:I][1]
    P[:t_s_c] = P[:t_s].-P[:t_s][I1]
    P[:F_kN_c] = P[:F_kN] .-P[:F_kN][I1]
    P[:U_mm_c] = ((P[:U1_mm].+P[:U2_mm]).-(P[:U1_mm][I1]+P[:U2_mm][I1]))./2
    P[:U_mm_fc] = P[:U_mm_c] .-(P[:F_kN_c]*exp_info[:K_mm_kN])
    P[:ε] = P[:U_mm_fc]./exp_info[:L_mm]
    P[:Jr] = JR(P, exp_info)
    P[:F_kN_j] = P[:F_kN_c] .-P[:Jr]
    P[:σ_MPa] = P[:F_kN_c]./(0.25e-6π*exp_info[:d_mm]^2) .*1e-3
    P[:σ_MPa_j] = P[:F_kN_j]./(0.25e-6π*exp_info[:d_mm]^2) .*1e-3
end

"""
    H_mod(P,ds::Int64)

    Calculate tangent moduli and yield stress using robust numerical differentiation

    #Arguments
    - 'P' : dictionary information pertaining to the 'Murrell'
"""
function Youngs_mod!(P, exp_info)
    ε1 = Array{Float64,1}(undef,1)
    σ1 = Array{Float64,1}(undef,1)
    t1 = Array{Float64,1}(undef,1)
    I1 = exp_info[:I][1]
    I2 = exp_info[:I][2]
    ε = P[:ε][I1:I2]
    σ = P[:σ_MPa][I1:I2]
    t = P[:t_s][I1:I2]
    # Youngs modulus calculation
    while true
        fig = figure()
        plot(ε,σ)
        xlim(0,0.005)
        ylim(0,300)
        g_data = ginput(2)
        close(fig)
        IE = (findfirst(ε.>g_data[1][1]), findfirst(ε.>g_data[2][1]))
        P[:EE] = linfit(ε[IE[1]:IE[2]], σ[IE[1]:IE[2]])
        P[:E_GPa] = P[:EE][1]/1e3
        fig = figure()
        plot(ε[IE[1]:IE[2]],σ[IE[1]:IE[2]],ε[IE[1]:IE[2]],M(ε[IE[1]:IE[2]],P[:EE].param))
        title("Youngs' modulus ="*string(P[:E_GPa])*"GPa")
        println("Fit ok? (y/n)")
        if readline() == "y"
            close(fig)
            break
        else
            close(fig)
        end
    end
end

function YS!(P,exp_info)
    I1 = exp_info[:I][1]
    I2 = exp_info[:I][2]
    ε1 = Array{Float64,1}(undef,1)
    σ1 = Array{Float64,1}(undef,1)
    t1 = Array{Float64,1}(undef,1)
    ε = P[:ε][I1:I2-1000]
    σ = P[:σ_MPa][I1:I2-1000]
    t = P[:t_s][I1:I2-1000]
    while true
        println("Enter moving average:")
        ds = parse(Float64,readline())
        for i in 1:Int(floor(length(P[:ε][I1:I2-1000])/ds))-1
            push!(ε1,sum(ε[Int(i*ds-ds+1):Int(i*ds)])/ds)
            push!(σ1,sum(σ[Int(i*ds-ds+1):Int(i*ds)])/ds)
            push!(t1,t[Int(i*ds)+Int(ceil(ds/2))])
        end
        deleteat!(ε1,1)
        deleteat!(σ1,1)
        deleteat!(t1,1)
        dσ = TVRegDiff(σ1,200, 0.2, dx=1, ε=1e-6)
        dε = TVRegDiff(ε1,200, 0.2, dx=1, ε=1e-6)
        P[:h] = dσ./dε*1e-3 # Calculate the tangent modulus in GPa, excluding initial non-linearity
        P[:ε_h] = ε1
        h_max = maximum(P[:h])
        II = findfirst(P[:h] .== h_max)
        P[:H_GPa] = sum(P[:h][P[:h].<0.4P[:E_GPa]])/sum(P[:h].<0.4P[:E_GPa])
        P[:ε_E] = ε1[II]
        P[:σ_E] = σ1[II]
        P[:YSa_MPa] = σ1[findfirst((P[:h].<0.8P[:E_GPa]) .& (ε1 .> ε1[II]))] #80% yield
        P[:εa] = ε1[findfirst((P[:h].<0.8P[:E_GPa]) .& (ε1 .> ε1[II]))]
        P[:YSb_MPa] = σ1[findfirst((P[:h].<0.5P[:E_GPa]) .& (ε1 .> ε1[II]))] #50% yield
        P[:εb] = ε1[findfirst((P[:h].<0.5P[:E_GPa]) .& (ε1 .> ε1[II]))]
        M(ε, p) = p[1]*ε .+p[2]
        fig = figure()
        ax=[]
        push!(ax,subplot(2,1,1))
        ax[1].plot(ε,σ,"0.5")
        ax[1].plot(ε,M(ε,P[:EE].param),"--k")
        ax[1].scatter(P[:εa],P[:YSa_MPa],s=100)
        ax[1].scatter(P[:εb],P[:YSb_MPa],s=100)
        ax[1].set_xlim(0,P[:εb]+0.001)
        ax[1].set_ylim(0,P[:YSb_MPa]+20)
        push!(ax,subplot(2,1,2))
        ax[2].plot(P[:ε_h],P[:h])
        ax[2].plot(P[:ε_h],M(P[:ε_h],(0,P[:E_GPa])),"--k")
        ax[2].plot(P[:ε_h],M(P[:ε_h],(0,0.8P[:E_GPa])),"--k")
        ax[2].plot(P[:ε_h],M(P[:ε_h],(0,0.5P[:E_GPa])),"--k")
        ax[2].scatter(P[:εa],0.8P[:E_GPa])
        ax[2].scatter(P[:εb],0.5P[:E_GPa])
        ax[2].set_xlim(0,P[:εb]+0.001)
        ax[2].set_ylim(0,2P[:E_GPa])
        println("Fit ok (y/n)")
        if readline() == "y"
            close(fig)
            break
        else
            close(fig)
        end
    end
end

"""
    SFP_read()

    Read and return data from program SFP_furnace_control_V2.vi

    #Arguments
    - 'SFP' : empty dictionary to store extracted data
    - 'fil::string' : path to file
"""
function SFP_read(fid)
    path="/Users/christopherharbord/Dropbox/My PC (DESKTOP-8JF2H49)/Documents/UCL/Furnace_calibration/SFP_logging/"*fil
    dat = readdlm(path)
    headers = dat[:,1]
    I = findall(x->x==headers[1],dat[:,1])
    SFP = Dict()
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
    JR(P, exp_info)

    Calculate strength contribution of copper jacket for corrections

    #Arguments
    * 'P' : dictionary containing processed data from the 'Murrell'
    * exp_info : details of experiment
"""
function JR!(P, exp_info)
    EaJ = 197000 # Copper activation enthalpy
    n1J = 4.8 # Empirical factor 1
    n2J = 0.22 # Empirical factor 2
    ε̇j = P[:εr]*(exp_info[:d_mm]/exp_info[:L_mm]) # Jacket strain rate
    Jr = 2*10^((log10(ε̇j*exp(EaJ/(8.3145*(exp_info[:T]+278)))))/n1J-n2J) # Copper flow stress at experiment conditions
    Ja = π*exp_info[:d_mm]*exp_info[:L_mm]*1e-6 # Jacket area
    JR = Jr.*Ja.*(P[:ε]*(exp_info[:d_mm]/exp_info[:L_mm])).*1e-3 # Force due to jacket assuming linear increase due to incremental strain
end

function M_interp!(P, t_UT)
    P[:F_kN_i] = lininterp(P[:t_s],P[:F_kN_j], t_UT)
    P[:U_mm_i] = lininterp(P[:t_s],P[:U_mm_c], t_UT)
    P[:σ_MPa_i] = lininterp(P[:t_s],P[:σ_MPa_j], t_UT)
    P[:σ3_MPa_i] = lininterp(P[:t_s],P[:PC2_MPa], t_UT)
    P[:ε_i] = lininterp(P[:t_s],P[:ε], t_UT)
end
