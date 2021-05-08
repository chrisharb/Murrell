module Murrell
"""
    A set of functions used to process data from the 'Murrell'
    gas medium triaxial, RIPL, UCL
"""
using TDMSReader, XLSX, NoiseRobustDifferentiation, Dates
using DelimitedFiles, OldTools

export M_file, M_read, M_reduce!, Youngs_mod!, YS!, SFP_read, JR!, M_interp!

include("base.jl")

end
