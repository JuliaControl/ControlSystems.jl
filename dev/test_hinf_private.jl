#@testset "test_hinf" begin

using ControlSystems
using Test
using LinearAlgebra

"""
Testing that we can re-write the weight specificaitons to a feasible state-space
form
"""
possible_entries = [[],1.0]

for ii = 1:2
  for kk = 1:2
    for ll = 1:2
      A   = rand(ii,ii)
      B   = rand(ii,kk)
      C   = rand(ll,ii)
      D   = rand(ll,kk)
      Gss = ss(A,B,C,D)
      Gtf = tf(Gss)
      append!(possible_entries, [Gss])
      append!(possible_entries, [Gtf])
    end
  end
end
#ii, jj, kk, ll= [4, 3, 2, 2]
for ii = 1:length(possible_entries)
  for jj = 2:length(possible_entries)
    for kk = 2:length(possible_entries)
      for ll = 2:length(possible_entries)
        G  = possible_entries[ii]
        WS = possible_entries[jj]
        WU = possible_entries[kk]
        WT = possible_entries[ll]
        println([ii, jj, kk, ll])
        # Check that the system is realizable
        if isa(G, StateSpace) || isa(G, TransferFunction) || isa(G, Number)
          Ag, Bg, Cg, Dg = ss(G).A, ss(G).B, ss(G).C, ss(G).D
          if (size(Ag,1)==0 || size(Ag,2)==0 || size(Bg,1)==0 || size(Bg,2)==0 ||
              size(Cg,1)==0 || size(Cg,2)==0 || size(Dg,1)==0 || size(Dg,2)==0)
            @test_throws ErrorException hInf_partition(G, WS, WU, WT)
          end
        elseif isempty(G)
          @test_throws ErrorException hInf_partition(G, WS, WU, WT)
        elseif isa(WT, StateSpace) || isa(WT, TransferFunction)
          if size(ss(G).C, 1) != size(ss(WT).B,2) && size(ss(WT).B,2) != 0
            @test_throws ErrorException hInf_partition(G, WS, WU, WT)
          end
        elseif isa(WS, StateSpace) || isa(WS, TransferFunction)
          if size(ss(G).A, 1) != size(ss(WS).B,2) && size(ss(WS).B,2) != 0
            @test_throws ErrorException hInf_partition(G, WS, WU, WT)
          end
        elseif isa(WU, StateSpace) || isa(WU, TransferFunction)
          if size(ss(G).B, 2) != size(ss(WU).B,2) && size(ss(WU).B,2) != 0
            @test_throws ErrorException hInf_partition(G, WS, WU, WT)
          end
        else
          A, B1, B2, C1, C2, D11, D12, D21, D22 = hInf_partition(G, WS, WU, WT)
          @test isa(ss(A,[B1 B2],[C1; C2], [D11 D12; D21 D22]), StateSpace)
        end
      end
    end
  end
end
