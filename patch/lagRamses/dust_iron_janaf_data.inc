! NIST-JANAF Fourth Edition (1998), pure Fe(cr,liquid), printed p1225.
! https://janaf.nist.gov/pdf/JANAF-FourthEd-1998-Iron.pdf (PDF page 7).
! Source SHA256: 96ec56da8d92685cf963a7382435064ef108d981cb36525b09c78fbf06213bbb
! Visually transcribed T[K], H(T)-H(298.15K)[kJ/mol]. Source mass 55.847 g/mol.
! Repeated T with different H is LATENT HEAT, not a malformed duplicate.
! The Curie point's identical H is stored once (no latent jump at 1042 K).
! Bulk equilibrium condensed phases only, not the FeO table or gas branch.
integer,parameter :: iron_n=40
real(real64),parameter :: iron_molar_mass=55.847d0,iron_h_zero=-4.507d0
real(real64),parameter :: iron_t(iron_n)=[ &
     0d0,100d0,200d0,298.15d0,300d0,400d0,500d0,600d0,700d0,800d0, &
     900d0,1000d0,1042d0,1100d0,1184d0,1184d0,1200d0,1300d0,1400d0,1500d0, &
     1600d0,1665d0,1665d0,1700d0,1800d0,1809d0,1809d0,1900d0,2000d0,2100d0, &
     2200d0,2300d0,2400d0,2500d0,2600d0,2700d0,2800d0,2900d0,3000d0,3100d0]
real(real64),parameter :: iron_h(iron_n)=[ &
     -4.507d0,-4.084d0,-2.315d0,0d0,.046d0,2.674d0,5.527d0,8.613d0,11.939d0,15.553d0, &
     19.582d0,24.329d0,26.960d0,29.970d0,33.619d0,34.519d0,35.062d0,38.505d0,42.033d0,45.643d0, &
     49.338d0,51.784d0,52.621d0,54.066d0,58.263d0,58.645d0,72.452d0,76.641d0,81.243d0,85.845d0, &
     90.448d0,95.050d0,99.653d0,104.255d0,108.857d0,113.460d0,118.062d0,122.665d0,127.267d0,131.869d0]
integer,parameter :: iron_transition_index(3)=[15,22,26]
real(real64),parameter :: iron_transition_t(3)=[1184d0,1665d0,1809d0]
