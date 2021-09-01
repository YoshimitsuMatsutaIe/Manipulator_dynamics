function (ˍ₋out, u)
    #= C:\Users\Yoshi\.julia\packages\SymbolicUtils\4n4Gz\src\code.jl:282 =#
    #= C:\Users\Yoshi\.julia\packages\SymbolicUtils\4n4Gz\src\code.jl:283 =#
    begin
        #= C:\Users\Yoshi\.julia\packages\Symbolics\CCDaf\src\build_function.jl:369 =#
        #= C:\Users\Yoshi\.julia\packages\SymbolicUtils\4n4Gz\src\code.jl:329 =# @inbounds begin
                #= C:\Users\Yoshi\.julia\packages\SymbolicUtils\4n4Gz\src\code.jl:325 =#
                ˍ₋out[1] = (+)((getindex)(u, 1), (*)(-1, (getindex)(u, 3)))
                ˍ₋out[2] = (+)((^)((getindex)(u, 1), 2), (*)(-1, (getindex)(u, 2)))
                ˍ₋out[3] = (+)((getindex)(u, 2), (getindex)(u, 3))
                #= C:\Users\Yoshi\.julia\packages\SymbolicUtils\4n4Gz\src\code.jl:327 =#
                nothing
            end
    end
end