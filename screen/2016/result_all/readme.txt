all sensing area
n edge, n node, n obs_odparis
"186;215;2343 "
car : 390
all object : 5523
markings : 1333
road sign : 188
obs kerb : 11690


odparis to odparis absolute distance :
mean median std
before optim
after optim 
1.856;1.507;1.534
0.481;0.356;0.469
0.281;0.118;0.442 -- no constraints


before/after kerb observation optimisation
odparis to edge surf absolute distance
1.856;1.507;1.534
0.825;0.608;0.783




before after object obs optim
distance between odparis and edge surf
1.856;1.507;1.534 
1.554;1.237;1.201


kerb and objects :
1.856;1.507;1.534 
0.878;0.656;0.791


--- all sensing area extended
only kerb : 

2.016;1.784;1.401
0.862;0.597;0.862

removing from measure edge with less than 30 observation each side
2.023;1.808;1.387
0.814;0.566;0.803




---------------------
--small area
----------------------
n edge, nnode, n obs odparis
"62;64;475 "
n obs kerb : 2628
user override : 62
cars : 145
objects : 2011 
signs : 344
marking : 136

optimising on kerb
distance between odparis and edge_surf_ext_ring before and after optim
2.044 & 1.8 & 0.937
0.986 & 0.533 & 1.331

same, but with added user override 
2.044 & 1.8 & 0.937
0.652 & 0.363 & 0.995

optimising on kerb and objects
2.044 & 1.8 & 0.937
1.026 & 0.716 & 1.117

kerb an dmakring
2.044 & 1.8 & 0.937
0.992 & 0.667 & 1.161

only marking
2.044 & 1.8 & 0.937
1.53 & 1.169 & 1.272

markings and cars
2.044 & 1.8 & 0.937
1.344 & 0.963 & 1.311


2.044 & 1.8 & 0.937
1.053 & 0.608 & 1.283

optimising using only odparis
2.044 & 1.8 & 0.937
0.677 & 0.431 & 1.055

----------------------------------------

for all paris:-- 4 min
1.821 & 1.565 & 1.379
0.722;0.526;0.759 
0.679;0.48;0.744 -- less constraints

all paris no constraint 
1.821 & 1.565 & 1.379
0.368 & 0.141 & 0.658


for whole paris :
edge node obs
"38821;45005;522270 "




------------------------------
very small area : rjmcmc markings
1.652;1.592;0.906
1.009;0.782;0.944










