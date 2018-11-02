using ControlSystems
using Test
using LinearAlgebra

@testset "test_hinf_bilinear" begin


  # Test SS-type data
  N = [1,5,10]
  M = [1,5,10]
  H = [0.1,0.01,0.001]

  for (ii, m) in enumerate(N)
    for (jj, n) in enumerate(M)
      for (kk, h) in enumerate(H)

        freq = [10^i for i in range(-6, stop=6, length=1001)]

        AcTrue, BcTrue, CcTrue, DcTrue = rand(m,m), rand(m,n), rand(n,m), rand(n,n)

        ev = eigvals(AcTrue)
        if !isempty(ev[real(ev).<0])
          AcTrue = AcTrue - one(AcTrue) * 2 * maximum(abs.(real(ev[real(ev).<0])))
        end

        Gtrue = ss(AcTrue, BcTrue, CcTrue, DcTrue)

        valA, fA = sigma(Gtrue, freq);
        sysB = hInf_bilinear_s2z(Gtrue, 0.05);
        valB, fB = sigma(sysB, freq)
        sysC = hInf_bilinear_z2s(hInf_bilinear_s2z(Gtrue, 0.05))
        valC, fC = sigma(sysC, freq)
        sysD = hInf_bilinear_s2z(hInf_bilinear_z2s(hInf_bilinear_s2z(Gtrue, 0.05)),0.05)
        valD, fD = sigma(sysD, freq)

        @test abs(maximum(svd(sysB.B).S)-maximum(svd(sysB.C).S)) < 1e-7
        @test abs(maximum(svd(sysC.B).S)-maximum(svd(sysC.C).S)) < 1e-7
        @test abs(maximum(svd(sysD.B).S)-maximum(svd(sysD.C).S)) < 1e-7

        # Test that the C->D->C and D->C->D results in the same
        @test norm(valA-valC) < 1e-7
        @test norm(valB-valD) < 1e-7
      end
    end
  end


  # Test extended ss-type data
  N  = [1,4]
  P1 = [1,4]
  P2 = [1,4]
  M1 = [1,4]
  M2 = [1,4]

  H = [0.1,0.01,0.001]

  for (in1, n) in enumerate(N)
    for (im1, m1) in enumerate(M1)
      for (im2, m2) in enumerate(M2)
        for (ip1, p1) in enumerate(P1)
          for (ip2, p2) in enumerate(P2)
            for (ih, h) in enumerate(H)

              freq = [10^i for i in range(-6, stop=6, length=1001)]

              A,  B1,  B2  = rand(n,n),  rand(n,m1),  rand(n,m2)
              C1, D11, D12 = rand(p1,n), rand(p1,m1), rand(p1,m2)
              C2, D21, D22 = rand(p2,n), rand(p2,m1), rand(p2,m2)

              ev = eigvals(A)
              if !isempty(ev[real(ev).>0])
                A = A - one(A) * 2 * maximum(abs.(real(ev[real(ev).>0])))
              end

              EsysA = ss(A, B1, B2, C1, C2, D11, D12, D21, D22)
              valA, fA = sigma(ss(
                  EsysA.A,
                  [EsysA.B1 EsysA.B2],
                  [EsysA.C1; EsysA.C2],
                  [EsysA.D11 EsysA.D12; EsysA.D21 EsysA.D22]) , freq)

              EsysB = hInf_bilinear_s2z(EsysA, h)
              valB, fB = sigma(ss(
                  EsysB.A,
                  [EsysB.B1 EsysB.B2],
                  [EsysB.C1; EsysB.C2],
                  [EsysB.D11 EsysB.D12; EsysB.D21 EsysB.D22]) , freq)

              EsysC = hInf_bilinear_z2s(hInf_bilinear_s2z(EsysA, h))
              valC, fC = sigma(ss(
                  EsysC.A,
                  [EsysC.B1 EsysC.B2],
                  [EsysC.C1; EsysC.C2],
                  [EsysC.D11 EsysC.D12; EsysC.D21 EsysC.D22]) , freq)

              EsysD = hInf_bilinear_s2z(hInf_bilinear_z2s(hInf_bilinear_s2z(EsysA, h)),h)
              valD, fD = sigma(ss(
                  EsysD.A,
                  [EsysD.B1 EsysD.B2],
                  [EsysD.C1; EsysD.C2],
                  [EsysD.D11 EsysD.D12; EsysD.D21 EsysD.D22]) , freq)

              # Test that the C->D->C and D->C->D results in the same
              @test norm(valA-valC) < 1e-7
              @test norm(valB-valD) < 1e-7
            end
          end
        end
      end
    end
  end
end

#plot(fA, valA, xscale = :log10, yscale = :log10, color = :blue, w=4)
#plot!(fB, valB, xscale = :log10, yscale = :log10, color = :red, w=3)
#plot!(fC, valC, xscale = :log10, yscale = :log10, color = :green, w=2)
#plot!(fD, valD, xscale = :log10, yscale = :log10, color = :black, w=1)
