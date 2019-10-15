mol addfile replace_prmtop
mol addfile mdcrdfile  first 0 last -1 step 1 waitfor -1
set nframe [molinfo top get numframes]
mol modstyle 0 0 NewCartoon 0.1 2
mol modcolor 0 0 Occupancy
mol modmaterial 0 0 Glass3
mol drawframes 0 0 0:$nframe
color Display Background white
display projection orthographic
display depthcue off
axes location off

mol addrep 0
animate goto need_to_decide
mol modstyle 1 0 NewCartoon
#mol modcolor 1 0 ColorID 1
mol modcolor 1 0 Structure
mol modmaterial 1 0 Occupancy

proc allign {{mol top}} {
    set reference [atomselect $mol "backbone" frame 0]
    set compare [atomselect $mol "backbone"]
    set num_steps [molinfo $mol get numframes]

    for {set frame 0} {$frame < $num_steps} {incr frame} {
    $compare frame $frame
    set trans_mat [measure fit $compare $reference]
    $compare move $trans_mat
    }
}

allign
