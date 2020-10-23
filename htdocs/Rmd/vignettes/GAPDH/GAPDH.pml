# @/Users/juancarlosaledo/ptm_outdropbox/ptm/htdocs/Rmd/vignettes/GAPDH/GAPDH.pml 


fetch 1u8f, async=0


# Esconder todo para luego ir mostrando de forma selectiva lo que queremos manejar.
hide everything, all
util.cbc # color por cadenas
set transparency, 0.5
set cartoon_side_chain_helper, on

## ------------- Cadena O ------------- ##

create obj_o, chain O
create MO46, (resi 46 and chain O)
create YO42, (i. 42 and chain O)
create YO94, (i. 94 and chain O)
create SO98, (i. 98 and chain O)
create TO99, (i. 99 and chain O)
create SO122, (i. 122 and chain O)
create TO237, (i. 237 and chain O)
create TO246, (i. 246 and chain O)

show ribbon, obj_o
color green, obj_o

show sticks, YO42
util.cbag YO42
show mesh, YO42

show sticks, MO46
util.cbag MO46
show mesh, MO46

show sticks, YO94
util.cbag YO94
show mesh, YO94

show sticks, SO98
util.cbag SO98
show mesh, SO98

show sticks, TO99
util.cbag TO99
show mesh, TO99

show sticks, SO122
util.cbag SO122
show mesh, SO122

show sticks, TO237
util.cbag TO237
show mesh, TO237

show sticks, TO246
util.cbag TO246
show mesh, TO246

create ptmO, (resi 42,46,94,98,99,122,237,246 and chain O)
show sticks, ptmO
util.cbag ptmO
show mesh, ptmO

## ------------- Cadena P ------------- ##
create obj_p, chain P
show ribbon, obj_p
color cyan, obj_p

create ptmP, (resi 42,46,94,98,99,122,237,246 and chain P)
show sticks, ptmP
util.cbac ptmP
show mesh, ptmP

## ------------- Cadena R ------------- ##

create obj_r, chain R
show ribbon, obj_r
color yellow, obj_r

create MR46, (resi 46 and chain R)
create YR42, (i. 42 and chain R)
create YR94, (i. 94 and chain R)
create SR98, (i. 98 and chain R)
create TR99, (i. 99 and chain R)
create SR122, (i. 122 and chain R)
create TR237, (i. 237 and chain R)
create TR246, (i. 246 and chain R)



create ptmR, (resi 42,46,94,98,99,122,237,246 and chain R)
show sticks, ptmR
util.cbay ptmR
show mesh, ptmR

# show sticks, MR46
# util.cbay MR46
# show mesh, MR46
# show sticks, YR42

# util.cbay YR42
# show mesh, YR42


