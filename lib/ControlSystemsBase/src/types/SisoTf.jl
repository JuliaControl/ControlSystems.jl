abstract type SisoTf{T<:Number} end

+(f1::SisoTf, f2::SisoTf) = +(promote(f1,f2)...)
-(f1::SisoTf, f2::SisoTf) = -(promote(f1,f2)...)

*(f1::SisoTf, f2::SisoTf) = *(promote(f1,f2)...)
/(f1::SisoTf, f2::SisoTf) = /(promote(f1,f2)...)
