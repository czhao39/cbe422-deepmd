LAMMPS (29 Oct 2020)
OMP_NUM_THREADS environment is not set. Defaulting to 1 thread. (src/comm.cpp:94)
  using 1 OpenMP thread(s) per MPI task
# Ar in metal units

# simulation params in reduced units
# settable from command line
# epsilon, sigma, mass set below

variable	x index 10
variable	y index 10
variable	z index 10
variable        rhostar index 0.8842
variable        dt index 0.005
variable        cutoff index 5.0
variable        skin index 0.3
variable        tinitial index 1.0
variable        nthermo index 10
variable        nsteps index 1000

# physical constants from update.cpp

variable        kb index 8.617343e-5          # kB in eV/K
variable        avogadro index 6.02214129e23  # Avogadro's number

# Ar properties in metal units

variable        epskb index 117.7             # LJ epsilon/kB in degrees K
variable        sigma index 3.504             # LJ sigma in Angstroms
variable        epsilon equal ${epskb}*${kb}  # LJ epsilon in eV
variable        epsilon equal 117.7*${kb}  
variable        epsilon equal 117.7*8.617343e-5  
variable        mass index 39.95              # mass in g/mole

# scale factors

# sigma = scale factor on distance, converts reduced distance to Angs
# epsilon = scale factor on energy, converts reduced energy to eV
# tmpscale = scale factor on temperature, converts reduced temp to degrees K
# tscale = scale factor on time, converts reduced time to ps
#   formula is t = t* / sqrt(epsilon/mass/sigma^2), but need t in fs
#   use epsilon (Joule), mass (kg/atom), sigma (meter) to get t in seconds
# pscale = scale factor on pressure, converts reduced pressure to bars
#   formula is P = P* / (sigma^3/epsilon), but need P in atmospheres
#   use sigma (meter), epsilon (Joule) to get P in nt/meter^2, convert to bars

variable        eVtoJoule index 1.602e-19     # convert eV to Joules
variable        NtMtoAtm equal 1.0e-5         # convert Nt/meter^2 to bars

variable        tmpscale equal ${epskb}
variable        tmpscale equal 117.7
variable        epsilonJ equal ${epsilon}*${eVtoJoule}
variable        epsilonJ equal 0.010142612711*${eVtoJoule}
variable        epsilonJ equal 0.010142612711*1.602e-19
variable        massKgAtom equal ${mass}/1000.0/${avogadro}
variable        massKgAtom equal 39.95/1000.0/${avogadro}
variable        massKgAtom equal 39.95/1000.0/6.02214129e23
variable        sigmaM equal ${sigma}/1.0e10
variable        sigmaM equal 3.504/1.0e10
variable        sigmaMsq equal ${sigmaM}*${sigmaM}
variable        sigmaMsq equal 3.504e-10*${sigmaM}
variable        sigmaMsq equal 3.504e-10*3.504e-10
variable        tscale equal 1.0e12/sqrt(${epsilonJ}/${massKgAtom}/${sigmaMsq})
variable        tscale equal 1.0e12/sqrt(1.6248465563022e-21/${massKgAtom}/${sigmaMsq})
variable        tscale equal 1.0e12/sqrt(1.6248465563022e-21/6.6338529895236e-26/${sigmaMsq})
variable        tscale equal 1.0e12/sqrt(1.6248465563022e-21/6.6338529895236e-26/1.2278016e-19)
variable        sigmaM3 equal ${sigmaM}*${sigmaM}*${sigmaM}
variable        sigmaM3 equal 3.504e-10*${sigmaM}*${sigmaM}
variable        sigmaM3 equal 3.504e-10*3.504e-10*${sigmaM}
variable        sigmaM3 equal 3.504e-10*3.504e-10*3.504e-10
variable        pscale equal ${NtMtoAtm}/(${sigmaM3}/(${epsilonJ}))
variable        pscale equal 1e-05/(${sigmaM3}/(${epsilonJ}))
variable        pscale equal 1e-05/(4.3022168064e-29/(${epsilonJ}))
variable        pscale equal 1e-05/(4.3022168064e-29/(1.6248465563022e-21))

# variables
# alat = lattice constant in Angs (at reduced density rhostar)
# temp = reduced temperature for output
# epair,emol,etotal = reduced epair,emol,etotal energies for output
# press = reduced pressure for output

variable        alat equal (4.0*${sigma}*${sigma}*${sigma}/${rhostar})^(1.0/3.0)
variable        alat equal (4.0*3.504*${sigma}*${sigma}/${rhostar})^(1.0/3.0)
variable        alat equal (4.0*3.504*3.504*${sigma}/${rhostar})^(1.0/3.0)
variable        alat equal (4.0*3.504*3.504*3.504/${rhostar})^(1.0/3.0)
variable        alat equal (4.0*3.504*3.504*3.504/0.8842)^(1.0/3.0)
variable        temp equal temp/${tmpscale}
variable        temp equal temp/117.7
variable        epair equal epair/${epsilon}
variable        epair equal epair/0.010142612711
variable        emol equal emol/${epsilon}
variable        emol equal emol/0.010142612711
variable        etotal equal etotal/${epsilon}
variable        etotal equal etotal/0.010142612711
variable        press equal press/${pscale}
variable        press equal press/377.676586146256

# same script as in.ar.lj

units		metal
atom_style	atomic

lattice		fcc ${alat}
lattice		fcc 5.79518437579763
Lattice spacing in x,y,z = 5.7951844 5.7951844 5.7951844
region		box block 0 $x 0 $y 0 $z
region		box block 0 10 0 $y 0 $z
region		box block 0 10 0 10 0 $z
region		box block 0 10 0 10 0 10
create_box	1 box
Created orthogonal box = (0.0000000 0.0000000 0.0000000) to (57.951844 57.951844 57.951844)
  1 by 2 by 2 MPI processor grid
create_atoms	1 box
Created 4000 atoms
  create_atoms CPU = 0.002 seconds
mass		1 ${mass}
mass		1 39.95

velocity	all create $(v_tinitial*v_epskb) 12345
velocity	all create 117.70000000000000284 12345

pair_style	deepmd ../shallow_metal.pb
pair_coeff

neighbor	$(v_skin*v_sigma) bin
neighbor	1.0511999999999999122 bin
neigh_modify	delay 0 every 20 check no

fix		1 all nve

timestep	$(v_dt*v_tscale)
timestep	0.011194658410003900315

# columns 2,3,4 = temp,pe,press in metal units
# columns 5-9 = temp,energy.press in reduced units, compare to in.ar.lj
# need to include metal unit output to enable use of reduced variables

thermo_style    custom step temp pe press v_temp v_epair v_emol v_etotal v_press
thermo_modify	norm yes
thermo		${nthermo}
thermo		10

run		${nsteps}
run		1000
Neighbor list info ...
  update every 20 steps, delay 0 steps, check no
  max neighbors/atom: 2000, page size: 100000
  master list distance cutoff = 18.5712
  ghost atom cutoff = 18.5712
  binsize = 9.2856, bins = 7 7 7
  1 neighbor lists, perpetual/occasional/extra = 1 0 0
  (1) pair deepmd, perpetual
      attributes: , newton on
      pair build: full/bin/atomonly
      stencil: full/bin/3d
      bin: standard
Per MPI rank memory allocation (min/avg/max) = 4.899 | 4.899 | 4.899 Mbytes
Step Temp PotEng Press v_temp v_epair v_emol v_etotal v_press 
       0        117.7  0.043538879    20606.034            1     4.292669            0     5.792294    54.559999 
      10    59.206406  0.051153789    21143.327   0.50302808    5.0434529            0    5.7978064    55.982626 
      20    64.471676  0.050569413    21048.324   0.54776275    4.9858369            0    5.8072757     55.73108 
      30    39.520949  0.053603659    21228.522   0.33577697    5.2849952            0    5.7885347    56.208201 
      40    57.221625  0.051331192    21084.761   0.48616504    5.0609437            0     5.790009    55.827557 
      50    49.685962  0.052268721    21145.946   0.42214071    5.1533784            0    5.7864312     55.98956 
      60    59.201991  0.051083584    21071.036   0.50299058    5.0365311            0    5.7908283    55.791216 
      70    57.654825  0.051277731    21089.892   0.48984558    5.0556728            0    5.7902575    55.841143 
      80    59.497656  0.051106989    21077.811    0.5055026    5.0388386            0    5.7969029    55.809155 
      90    58.989893  0.051192735    21089.228   0.50118856    5.0472927            0    5.7988876    55.839384 
     100    60.972235  0.050956189    21078.105   0.51803088    5.0239707            0    5.8008227    55.809932 
     110    61.250923  0.050966722    21077.033   0.52039867    5.0250092            0     5.805412    55.807096 
     120    61.979324  0.050881288    21072.586   0.52658729    5.0165859            0    5.8062693     55.79532 
     130    59.846578  0.051154338    21091.087    0.5084671     5.043507            0     5.806017    55.844307 
     140    61.044874  0.051007213    21079.193   0.51864803    5.0290014            0    5.8067789    55.812813 
     150    58.634352  0.051334607    21095.075   0.49816781    5.0612804            0    5.8083453    55.854867 
     160    60.092368  0.051149972    21087.267   0.51055538    5.0430765            0    5.8087181    55.834192 
     170      59.2405  0.051286787    21092.716   0.50331776    5.0565656            0    5.8113535    55.848619 
     180    59.592731  0.051253233    21089.877   0.50631037    5.0532574            0    5.8125331    55.841104 
     190    59.551623  0.051240972    21092.224   0.50596111    5.0520486            0    5.8108005    55.847318 
     200    60.531188  0.051126729    21086.518   0.51428367    5.0407849            0    5.8120176    55.832208 
     210    60.526713  0.051155336    21089.057   0.51424565    5.0436054            0     5.814781    55.838931 
     220    60.071822  0.051226117    21091.735   0.51038081    5.0505839            0    5.8159638    55.846021 
     230    61.578824  0.051062372    21083.056   0.52318457    5.0344397            0    5.8190204    55.823041 
     240    59.926166  0.051273522    21100.873    0.5091433    5.0552578            0    5.8187818    55.870219 
     250    61.449119  0.051100373    21087.972   0.52208257    5.0381863            0    5.8211144    55.836059 
     260    60.323449  0.051241609    21095.077   0.51251869    5.0521114            0    5.8206972    55.854871 
     270    61.475334  0.051106025    21083.869    0.5223053    5.0387436            0    5.8220057    55.825195 
     280    60.492283  0.051240302    21095.069   0.51395312    5.0519825            0    5.8227194    55.854851 
     290     60.68519  0.051231314    21094.339    0.5155921    5.0510964            0    5.8242912    55.852918 
     300    60.747402  0.051253291    21095.536   0.51612066    5.0532632            0    5.8272506    55.856086 
     310    61.252619  0.051190165    21092.309   0.52041308    5.0470393            0    5.8274637    55.847541 
     320    61.296147  0.051188145     21093.81    0.5207829    5.0468401            0    5.8278192    55.851515 
     330     60.61084  0.051283718    21101.463   0.51496041    5.0562631            0    5.8285106    55.871781 
     340    60.812822  0.051272752    21103.229   0.51667648    5.0551819            0    5.8300029    55.876455 
     350    61.278503  0.051244204    21101.878   0.52063299    5.0523672            0    5.8331214    55.872879 
     360    61.467458  0.051217954     21098.58   0.52223839    5.0497792            0    5.8329409    55.864147 
     370    61.852747  0.051196873    21096.835   0.52551187    5.0477006            0    5.8357714    55.859527 
     380    62.303329  0.051128645    21091.697    0.5293401    5.0409738            0    5.8347854    55.845922 
     390    61.931131  0.051221804    21097.793   0.52617784    5.0501587            0    5.8392282    55.862063 
     400    61.267117  0.051319538     21102.61   0.52053626    5.0597946            0    5.8404038    55.874818 
     410    60.644308  0.051389819    21108.345   0.51524476     5.066724            0    5.8393979    55.890002 
     420    61.466157  0.051288922    21102.878   0.52222733    5.0567761            0    5.8399213    55.875526 
     430    61.107989  0.051341078    21105.174   0.51918427    5.0619184            0    5.8405001    55.881606 
     440    61.074947  0.051362821     21106.86   0.51890354    5.0640621            0    5.8422228    55.886069 
     450    61.046535   0.05135742    21105.621   0.51866215    5.0635297            0    5.8413284     55.88279 
     460    61.909369  0.051250323    21098.899   0.52599294    5.0529705            0    5.8417627    55.864992 
     470    61.836011  0.051275784     21103.85   0.52536967    5.0554808            0    5.8433383    55.878099 
     480    62.484728  0.051211962    21100.358   0.53088129    5.0491884            0    5.8453112    55.868853 
     490    60.514251  0.051483573     21112.89   0.51413977    5.0759676            0    5.8469844    55.902035 
     500    61.643466  0.051340781    21105.251   0.52373378    5.0618891            0    5.8472934    55.881811 
     510     61.41651  0.051380645    21107.615   0.52180552    5.0658195            0    5.8483321    55.888069 
     520    61.724086  0.051348158    21105.551   0.52441874    5.0626165            0    5.8490479    55.882604 
     530    62.460061  0.051245983    21100.335   0.53067172    5.0525426            0    5.8483512    55.868794 
     540    61.170012  0.051422632     21111.46   0.51971123    5.0699592            0    5.8493311    55.898249 
     550    62.120465  0.051326544    21105.228   0.52778645    5.0604854            0    5.8519672    55.881748 
     560    61.441621  0.051410765    21109.523   0.52201887    5.0687891            0    5.8516217     55.89312 
     570    61.368126  0.051428122    21111.735   0.52139444    5.0705004            0    5.8523965    55.898978 
     580    62.538796  0.051313125     21104.67   0.53134066    5.0591624            0    5.8559741    55.880272 
     590    63.367849  0.051226166    21101.467   0.53838444    5.0505888            0    5.8579635    55.871791 
     600    63.690182  0.051183161    21098.348   0.54112304    5.0463487            0    5.8578304    55.863531 
     610    61.534618  0.051456739    21116.028   0.52280899    5.0733219            0    5.8573393    55.910344 
     620    61.549256  0.051470858    21119.299   0.52293335    5.0747139            0    5.8589179    55.919004 
     630    61.959407  0.051447887    21113.677   0.52641807    5.0724491            0    5.8618788    55.904119 
     640    62.862869  0.051341095    21109.524   0.53409404    5.0619201            0    5.8628609    55.893124 
     650    61.833817  0.051481889    21118.097   0.52535104    5.0758015            0     5.863631    55.915823 
     660    62.857597  0.051362995    21111.425   0.53404925    5.0640793            0    5.8649529    55.898157 
     670    61.886806  0.051513714    21119.934   0.52580124    5.0789393            0    5.8674439    55.920686 
     680    63.071585  0.051387577    21108.251   0.53586733    5.0665029            0    5.8701029    55.889753 
     690    62.288596     0.051485    21117.044   0.52921492    5.0761082            0    5.8697322    55.913036 
     700    63.806512  0.051300722    21107.657    0.5421114    5.0579395            0    5.8709033    55.888181 
     710    62.648574  0.051448312    21119.343   0.53227335     5.072491            0    5.8707015    55.919122 
     720    62.767185   0.05145791    21118.882    0.5332811    5.0734373            0     5.873159    55.917901 
     730    62.806359  0.051463618     21118.74   0.53361393    5.0740001            0    5.8742209    55.917525 
     740    61.501525  0.051621346    21127.099   0.52252782    5.0895511            0    5.8731469    55.939657 
     750    63.242359  0.051442672    21115.817   0.53731826    5.0719349            0    5.8777108    55.909787 
     760    62.627291  0.051516629    21119.492   0.53209253    5.0792267            0     5.877166    55.919517 
     770    63.651491   0.05141639    21114.202   0.54079432    5.0693437            0    5.8803324    55.905509 
     780    61.703164  0.051655166    21132.162   0.52424098    5.0928856            0    5.8790504    55.953063 
     790    62.824986  0.051518779    21122.804   0.53377218    5.0794386            0    5.8798967    55.928287 
     800    62.609432  0.051572825    21123.879    0.5319408    5.0847673            0     5.882479    55.931131 
     810    63.841108  0.051427216    21116.705   0.54240534    5.0704111            0    5.8838157    55.912137 
     820    63.702511  0.051441018    21120.803    0.5412278    5.0717718            0    5.8834106    55.922987 
     830    63.913463  0.051428242    21120.616   0.54302008    5.0705123            0    5.8848388    55.922493 
     840    63.841407  0.051462487    21122.284   0.54240788    5.0738886            0     5.887297    55.926908 
     850    63.330689  0.051527799     21124.17   0.53806873     5.080328            0    5.8872293    55.931903 
     860     63.73191  0.051528574    21122.941   0.54147757    5.0804044            0    5.8924177    55.928648 
     870    63.660862  0.051534517    21124.406   0.54087393    5.0809903            0    5.8920984    55.932527 
     880    62.697814  0.051647129    21131.702   0.53269171    5.0920932            0     5.890931    55.951845 
     890    64.314933  0.051447977    21120.623   0.54643104     5.072458            0    5.8918997     55.92251 
     900    63.036004   0.05162444    21133.396   0.53556503    5.0898562            0    5.8930029     55.95633 
     910    63.185336  0.051618071    21131.828   0.53683378    5.0892282            0    5.8942776    55.952178 
     920    63.749123  0.051545624    21128.842   0.54162381    5.0820854            0     5.894318    55.944273 
     930    64.580205  0.051465671    21120.703   0.54868483    5.0742025            0     5.897024    55.922723 
     940    63.425276  0.051624064    21131.032   0.53887236    5.0898191            0    5.8979256    55.950071 
     950    63.905054  0.051565988    21126.766   0.54294863    5.0840931            0    5.8983125    55.938776 
     960    63.520751  0.051651038    21131.143   0.53968352    5.0924786            0    5.9018015    55.950366 
     970     63.14566  0.051693056    21136.318   0.53649669    5.0966213            0    5.9011651    55.964068 
     980    63.762809  0.051605583     21132.89   0.54174009     5.087997            0     5.900404     55.95499 
     990     64.20021  0.051565524    21131.513   0.54545633    5.0840474            0    5.9020273    55.951346 
    1000    64.964322  0.051497651     21126.31   0.55194836    5.0773556            0    5.9050711    55.937568 
Loop time of 255.047 on 4 procs for 1000 steps with 4000 atoms

Performance: 3.792 ns/day, 6.329 hours/ns, 3.921 timesteps/s
186.2% CPU use with 4 MPI tasks x 1 OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 239.86     | 240.96     | 241.99     |   5.6 | 94.48
Neigh   | 2.3715     | 2.4348     | 2.4839     |   3.0 |  0.95
Comm    | 10.488     | 11.55      | 12.685     |  25.5 |  4.53
Output  | 0.012183   | 0.016001   | 0.018849   |   2.2 |  0.01
Modify  | 0.062234   | 0.062731   | 0.063248   |   0.1 |  0.02
Other   |            | 0.02019    |            |       |  0.01

Nlocal:        1000.00 ave        1009 max         994 min
Histogram: 2 0 0 0 0 0 1 0 0 1
Nghost:        7729.00 ave        7735 max        7720 min
Histogram: 1 0 0 0 1 0 0 0 0 2
Neighs:        0.00000 ave           0 max           0 min
Histogram: 4 0 0 0 0 0 0 0 0 0

Total # of neighbors = 0
Ave neighs/atom = 0.0000000
Neighbor list builds = 50
Dangerous builds not checked

Please see the log.cite file for references relevant to this simulation

Total wall time: 0:04:15
