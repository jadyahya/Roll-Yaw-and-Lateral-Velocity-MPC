syms m vdot u r ms hs phidd Fy ixx phi g ks ls bs phid phir
eqn1 = m*(vdot+u*r-g*phir)+ms*hs*phidd == Fy;
syms izz rd Mz
eqn2 = izz*rd == Mz;
eqn3 = (ixx+m*hs^2)*phidd == -m*(vdot+u*r-g*phir)*hs+m*g*hs*(phi-phir)-0.5*ks*ls^2*(phi)-0.5*bs*ls^2*phid;
S = solve([eqn1 eqn2 eqn3],[vdot rd phidd])