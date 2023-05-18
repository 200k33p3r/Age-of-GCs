#! /bin/bash
#

setups ()
{


#set up names of various files and standard parmeters for newpoly

# this files will not change
    ln -s $opac/FERMI.TAB                  fort.15

#setup a bunch of standard parameters for the polytrope program
    beta=1.000
    ddage=10000.0
    fmass1=0.0000099
    fmass2=0.9999
    pn=1.5
    lexcom=.TRUE.
}


############################################################

get_teff_lum ()
{
   if (( $mass <= 390 )); then
        teff=$(awk "BEGIN {print 0.430 * $sumass + 3.410}")
        lum=$(awk "BEGIN {print  6.20 * $sumass  -1.20}")
    else
#	teff=3.64
#	teff=$(awk "BEGIN {print 3.64 +  (1.0 - $sumass)*0.05}")
#	lum=$(awk "BEGIN {print 2.6 + 0.2 * ($sumass - 5.0)}")
       teff=3.88
       teff=$(awk "BEGIN {print 3.88 + (1.0 - $sumass)*0.05}")
       lum=$(awk "BEGIN {print 2.9 + 0.2 * ($sumass-5.0)}")
    fi
#    teff=$(awk "BEGIN {print $teff + $Tadjust}")
#    lum=$(awk "BEGIN {print $lum + $Ladjust}")
    teff=$(awk "BEGIN {print $Tadjust}")
    lum=$(awk "BEGIN {print $Ladjust}")

}

############################################################

wrtpoly ()
{

#creat the name list file for newpoly.f
    if [[ -e poly.nml ]]; then
	rm poly.nml
    fi
    echo "\$DATA " >poly.nml
    echo "SUMASS= $sumass" >>poly.nml
    echo "TEFFL1= $teff" >>poly.nml
    echo "SULUML= $lum" >>poly.nml
    echo "X= $xx" >>poly.nml
    echo "Z= $zz" >>poly.nml
    echo "ELEM(1)=${elem[1]}" >>poly.nml
    for ((j=2; j <= 8; j++ )); do
	z_elem=$(awk "BEGIN {print $zz*${elem[$j]}}")
	echo "ELEM($j)=$z_elem" >>poly.nml
    done
    echo "CMIXL= $mixlen" >>poly.nml
    echo "BETA= $beta" >>poly.nml
    echo "FMASS1= $fmass1" >>poly.nml
    echo "FMASS2= $fmass2" >>poly.nml
    echo "DDAGE= $ddage" >>poly.nml
    echo "PN= $pn" >>poly.nml
    echo "LEXCOM=$lexcom" >>poly.nml
    echo "\$END" >>poly.nml

}
	
############################################################

getphysnml ()
{
    if [[ -e physics.nml ]]; then
	rm physics.nml
    fi
    echo "\$PHYSICS " >physics.nml
    if [[ $mass -lt 691 ]]; then
	cat $nml/phys1.low.nml >>physics.nml
    elif [[ $mass -lt 1501 ]]; then
	cat $nml/phys1medm.nml >>physics.nml
    else
        cat $nml/physhighm.nml >>physics.nml
    fi
    ln -s physics.nml fort.13
}

############################################################

setupmix ()
{

#set up the namelist, the opacity and the surface boundary conditions
#files for the specific mixture (mix ID, [Fe/H] and eventually [a/Fe] )
# defaults are based upon NACRE rates & NO diffussion; must change
# zsolar, mixlen, xx (solar) for other physics

    if [[ -e control.nml ]]; then
	rm control.nml
    fi
	
    if [[ $ID = "GAS07" ]] ; then
	if [[ "$afe" != "0.0" ]]; then
	    echo "You specified [a/Fe] = $afe"
	    echo "for GAS07 mixture, can only do [a/Fe] = 0.0"
	    exit 1
	fi
	elem=(1.0E22 3.0E-5 0.17533 0.00177 0.050696 0.0 0.439279 0.0 0.0)
	zsolar=0.01405
#calculate z from [Fe/H]
	zz=$(awk "BEGIN {print (10.0^$feh)*$zsolar}")

	echo "zz= $zz   [Fe/H] = $feh"
#setup switch for boundry conditions
	bcswitch=noalpha

#	need to creat the name list file
	if [[ $mixlen = "solar" ]]; then
	    mixlen=1.9126
	fi
	if [[ $xx = "solar" ]]; then
	    xx=0.72520
	fi
	cat $nml/gas07_1.nml >control.nml
	echo " NMODLS(2)=$nmod" >>control.nml
	echo " RSCLZ(1)=$zz " >>control.nml
	echo " RSCLX(1)=$xx " >>control.nml
	echo " CMIXLA(1)=$mixlen " >>control.nml
	echo " CMIXLA(2)=$mixlen " >>control.nml
	echo " ENDAGE(2)=$endage " >>control.nml
	echo "  " >>control.nml
	echo " ZALEX=$zz" >>control.nml
	cat $nml/gas07_2.nml >>control.nml

	ln -s control.nml                    fort.14
	ln -s $opal/AGS04hz                  fort.48

    elif [[ $ID = "GS98" ]] ; then 
	zsolar=0.01889
	if [[ $mixlen = "solar" ]]; then
	    mixlen=1.9410
	fi
	if [[ $xx = "solar" ]]; then
	    xx=0.70692
	fi

# create first part of the control name list file
	cat $nml/gs98_1.nml >control.nml
	echo " NMODLS(2)=$nmod" >>control.nml
	echo " RSCLZ(1)=$zz " >>control.nml
	echo " RSCLX(1)=$xx " >>control.nml
	echo " CMIXLA(1)=$mixlen " >>control.nml
	echo " CMIXLA(2)=$mixlen " >>control.nml
	echo " ENDAGE(2)=$endage " >>control.nml
	echo "  " >>control.nml
	echo " ZALEX=$zz" >>control.nml
	echo "  " >>control.nml


	if [[ $afe = "0.000000" ]]; then
	    elem=(1.0E22 3.0E-5 0.17208 0.00150 0.050410 0.0 0.468020 0.0 0.0)
	    echo "in setupmix, we have picked [a/Fe]=0.0"
#setup switch for boundry conditions
	    bcswitch=noalpha

	    cat $nml/gs98_3_DATA.nml >>control.nml
	    
	    ln -s control.nml                    fort.14
	    ln -s $opal/GS98hz                   fort.48

	elif [[ $afe = "-0.200000" ]]; then
	    elem=(1.0E22 3.0E-5 0.22961 0.00150 0.067260 0.0 0.394030 0.0 0.0) 
	    echo "in setupmix, we have picked [a/Fe]=-0.2"

	    bcswitch=alphastand
	    bcsuf=m0d2

	    cat $nml/gs98afem2_3_DATA.nml >>control.nml
	    
	    ln -s control.nml                    fort.14
	    ln -s $opal/GS98hz_OFe-.2            fort.48
	elif [[ $afe = "0.200000" ]]; then
	    elem=(1.0E22 3.0E-5 0.12329 1.5E-3 0.03611 0.0 0.53143 0.0 0.0)
	    echo "in setupmix, we have picked [a/Fe]=+0.2"

	    bcswitch=alphastand
	    bcsuf=p0d2

	    cat $nml/gs98afep2_3_DATA.nml >>control.nml
	    
	    ln -s control.nml                    fort.14
	    ln -s $opal/GS98hz_OFe.2            fort.48
	elif [[ $afe = "0.400000" ]]; then
	    elem=(1.0E22 3.0E-5 0.08490 1.5E-3 0.02487 0.0 0.58003 0.0 0.0)
	    echo "in setupmix, we have picked [a/Fe]=+0.4"

	    bcswitch=alphastand
	    bcsuf=p0d4

	    cat $nml/gs98afep4_3_DATA.nml >>control.nml
	    
	    ln -s control.nml                    fort.14
	    ln -s $opal/GS98hz_OFe.4            fort.48
	elif [[ $afe = "0.600000" ]]; then
	    elem=(1.0E22 3.0E-5 0.05685 1.5E-3 0.01665 0.0 0.61554 0.0 0.0)
	    echo "in setupmix, we have picked [a/Fe]=+0.6"

	    bcswitch=alphastand
	    bcsuf=p0d6

	    cat $nml/gs98afep6_3_DATA.nml >>control.nml
	    
	    ln -s control.nml                    fort.14
	    ln -s $opal/GS98hz_OFe.6            fort.48
	elif [[ $afe = "0.800000" ]]; then
	    elem=(1.0E22 3.0E-5 0.03731 1.5E-3 0.01093 0.0 0.64027 0.0 0.0)
	    echo "in setupmix, we have picked [a/Fe]=+0.8"

	    bcswitch=alphastand
	    bcsuf=p0d8

	    cat $nml/gs98afep8_3_DATA.nml >>control.nml
	    
	    ln -s control.nml                    fort.14
	    ln -s $opal/GS98hz_OFe.8            fort.48
	else
	    echo "You specified [a/Fe] = $afe"
	    echo "For GS98, [a/Fe] = -0.2; 0.0, +0.2; +0.4; +0.6; +0.8 only"
	    exit 1
	fi

    else
	echo "Incorrect mix ID; specified $ID  only 'GAS07' or 'GS98' allowed"
	exit 1
    fi


#select the appropriate boundry condition file
#need to pick out closest [Fe/H] to the tabulated values
    fehint=$(awk "BEGIN {printf "'"%3d"'", $feh*100}")

    if [[ $bcswitch = "noalpha" ]]; then
	if (( $fehint <= -375 )); then
	    echo "[Fe/h] = -4.0 BC used"
	    ln -s $bckur/atmk1990m40.tab           fort.38
	    ln -s $bcphx/z_m4d0.afe_p0d0.dat       fort.95
	    ln -s $bcphx/z_m3d5.afe_p0d0.dat       fort.96
	    ln -s $bcphx/z_m3d0.afe_p0d0.dat       fort.97
	    ln -s $bcphx/z_m2d5.afe_p0d0.dat       fort.98
	    ln -s $bcphx/z_m2d0.afe_p0d0.dat       fort.99
	elif (( $fehint <= -325 )); then
	    echo "[Fe/H] = -3.5 BC used"
	    ln -s $bckur/atmk1990m35.tab           fort.38
	    ln -s $bcphx/z_m3d5.afe_p0d0.dat       fort.95
	    ln -s $bcphx/z_m3d0.afe_p0d0.dat       fort.96
	    ln -s $bcphx/z_m2d5.afe_p0d0.dat       fort.97
	    ln -s $bcphx/z_m2d0.afe_p0d0.dat       fort.98
	    ln -s $bcphx/z_m1d5.afe_p0d0.dat       fort.99
	elif (( $fehint <= -275 )); then
	    echo "[Fe/H] = -3.0 BC used"
	    ln -s $bckur/atmk1990m30.tab           fort.38
	    ln -s $bcphx/z_m3d0.afe_p0d0.dat       fort.95
	    ln -s $bcphx/z_m2d5.afe_p0d0.dat       fort.96
	    ln -s $bcphx/z_m2d0.afe_p0d0.dat       fort.97
	    ln -s $bcphx/z_m1d5.afe_p0d0.dat       fort.98
	    ln -s $bcphx/z_m1d0.afe_p0d0.dat       fort.99
	elif (( $fehint <= -225 )); then
	    echo "[Fe/H] = -2.5 BC used"
	    ln -s $bckur/atmk1990m25.tab           fort.38
	    ln -s $bcphx/z_m2d5.afe_p0d0.dat       fort.95
	    ln -s $bcphx/z_m2d0.afe_p0d0.dat       fort.96
	    ln -s $bcphx/z_m1d5.afe_p0d0.dat       fort.97
	    ln -s $bcphx/z_m1d0.afe_p0d0.dat       fort.98
	    ln -s $bcphx/z_m0d7.afe_p0d0.dat       fort.99
	elif (( $fehint <= -175 )); then
	    echo "[Fe/H] = -2.0 BC used"
	    ln -s $bckur/atmk1990m20.tab           fort.38
	    ln -s $bcphx/z_m2d0.afe_p0d0.dat       fort.95
	    ln -s $bcphx/z_m1d5.afe_p0d0.dat       fort.96
	    ln -s $bcphx/z_m1d0.afe_p0d0.dat       fort.97
	    ln -s $bcphx/z_m0d7.afe_p0d0.dat       fort.98
	    ln -s $bcphx/z_m0d5.afe_p0d0.dat       fort.99
	elif (( $fehint <= -125 )); then
	    echo "[Fe/H] = -1.5 BC used"
	    ln -s $bckur/atmk1990m15.tab           fort.38 
	    ln -s $bcphx/z_m1d5.afe_p0d0.dat       fort.95 
	    ln -s $bcphx/z_m1d0.afe_p0d0.dat       fort.96 
	    ln -s $bcphx/z_m0d7.afe_p0d0.dat       fort.97 
	    ln -s $bcphx/z_m0d5.afe_p0d0.dat       fort.98 
	    ln -s $bcphx/z_p0d0.afe_p0d0.dat       fort.99
	elif (( $fehint <= -85 )); then
	    echo "[Fe/H] = -1.0 BC used"
	    ln -s $bckur/atmk1990m10.tab           fort.38 
	    ln -s $bcphx/z_m1d0.afe_p0d0.dat       fort.95 
	    ln -s $bcphx/z_m0d7.afe_p0d0.dat       fort.96 
	    ln -s $bcphx/z_m0d5.afe_p0d0.dat       fort.97 
	    ln -s $bcphx/z_p0d0.afe_p0d0.dat       fort.98 
	    ln -s $bcphx/z_p0d15.afe_p0d0.dat      fort.99
	elif (( $fehint <= -60 )); then
	    echo "[Fe/H] = -0.7 BC used"
	    ln -s $bckur/atmk1990m07.tab           fort.38 
	    ln -s $bcphx/z_m0d7.afe_p0d0.dat       fort.95 
	    ln -s $bcphx/z_m0d5.afe_p0d0.dat       fort.96 
	    ln -s $bcphx/z_p0d0.afe_p0d0.dat       fort.97 
	    ln -s $bcphx/z_p0d15.afe_p0d0.dat      fort.98 
	    ln -s $bcphx/z_p0d3.afe_p0d0.dat       fort.99
	elif (( $fehint <= -25 )); then
	    echo "[Fe/H] = -0.5 BC used"
	    ln -s $bckur/atmk1990m05.tab           fort.38 
	    ln -s $bcphx/z_m0d5.afe_p0d0.dat       fort.95 
	    ln -s $bcphx/z_p0d0.afe_p0d0.dat       fort.96 
	    ln -s $bcphx/z_p0d15.afe_p0d0.dat      fort.97 
	    ln -s $bcphx/z_p0d3.afe_p0d0.dat       fort.98 
	    ln -s $bcphx/z_p0d5.afe_p0d0.dat       fort.99
	elif (( $fehint <= 8 )); then
	    echo "[Fe/H] = 0.0 BC used"
	    ln -s $bckur/atmk1990p00.tab           fort.38 
	    ln -s $bcphx/z_p0d0.afe_p0d0.dat       fort.95 
	    ln -s $bcphx/z_m0d5.afe_p0d0.dat       fort.96 
	    ln -s $bcphx/z_m0d7.afe_p0d0.dat       fort.97 
	    ln -s $bcphx/z_m1d0.afe_p0d0.dat       fort.98 
	    ln -s $bcphx/z_m1d5.afe_p0d0.dat       fort.99
	elif (( $fehint <= 22 )); then
	    echo "[Fe/H] = +0.15 BC used"
	    ln -s $bckur/atmk1990p015.tab          fort.38 
	    ln -s $bcphx/z_p0d15.afe_p0d0.dat      fort.95 
	    ln -s $bcphx/z_p0d0.afe_p0d0.dat       fort.96 
	    ln -s $bcphx/z_m0d5.afe_p0d0.dat       fort.97 
	    ln -s $bcphx/z_m0d7.afe_p0d0.dat       fort.98 
	    ln -s $bcphx/z_m1d0.afe_p0d0.dat       fort.99
	elif (( $fehint <= 40 )); then
	    echo "[Fe/H] = +0.3 BC used"
	    ln -s $bckur/atmk1990p03.tab           fort.38 
	    ln -s $bcphx/z_p0d3.afe_p0d0.dat       fort.95 
	    ln -s $bcphx/z_p0d15.afe_p0d0.dat      fort.96 
	    ln -s $bcphx/z_p0d0.afe_p0d0.dat       fort.97 
	    ln -s $bcphx/z_m0d5.afe_p0d0.dat       fort.98 
	    ln -s $bcphx/z_m0d7.afe_p0d0.dat       fort.99
	else
	    echo "[Fe/H] = +0.5 BC used"
	    ln -s $bckur/atmk1990p05.tab           fort.38 
	    ln -s $bcphx/z_p0d5.afe_p0d0.dat       fort.95 
	    ln -s $bcphx/z_p0d3.afe_p0d0.dat       fort.96 
	    ln -s $bcphx/z_p0d15.afe_p0d0.dat      fort.97 
	    ln -s $bcphx/z_p0d0.afe_p0d0.dat       fort.98 
	    ln -s $bcphx/z_m0d5.afe_p0d0.dat       fort.99
	fi
    elif [[ $bcswitch = "alphastand" ]]; then
#standard case for alpha-enhanced boundary condition files
	if (( $fehint <= -325 )); then
	    echo "[Fe/H] = -3.0 BC [a/Fe]=0.0 used"
	    ln -s $bckur/atmk1990m30.tab           fort.38
	    ln -s $bcphx/z_m3d0.afe_p0d0.dat       fort.95
	    ln -s $bcphx/z_m2d5.afe_p0d0.dat       fort.96
	    ln -s $bcphx/z_m2d0.afe_p0d0.dat       fort.97
	    ln -s $bcphx/z_m1d5.afe_p0d0.dat       fort.98
	    ln -s $bcphx/z_m1d0.afe_p0d0.dat       fort.99
	elif (( $fehint <= -225 )); then
	    echo "[Fe/H] = -2.5 BC used"
	    ln -s $bckur/atmk1990m25.tab           fort.38
	    ln -s $bcphx/z_m2d5.afe_$bcsuf.dat     fort.95
	    ln -s $bcphx/z_m2d0.afe_$bcsuf.dat     fort.96
	    ln -s $bcphx/z_m1d5.afe_$bcsuf.dat     fort.97
	    ln -s $bcphx/z_m1d0.afe_$bcsuf.dat     fort.98
	    ln -s $bcphx/z_m0d5.afe_$bcsuf.dat     fort.99
	elif (( $fehint <= -175 )); then
	    echo "[Fe/H] = -2.0 BC used"
	    ln -s $bckur/atmk1990m20.tab           fort.38
	    ln -s $bcphx/z_m2d0.afe_$bcsuf.dat     fort.95
	    ln -s $bcphx/z_m1d5.afe_$bcsuf.dat     fort.96
	    ln -s $bcphx/z_m1d0.afe_$bcsuf.dat     fort.97
	    ln -s $bcphx/z_m0d5.afe_$bcsuf.dat     fort.98
	    ln -s $bcphx/z_p0d0.afe_$bcsuf.dat     fort.99
	elif (( $fehint <= -125 )); then
	    echo "[Fe/H] = -1.5 BC used"
	    ln -s $bckur/atmk1990m15.tab           fort.38 
	    ln -s $bcphx/z_m1d5.afe_$bcsuf.dat     fort.95 
	    ln -s $bcphx/z_m1d0.afe_$bcsuf.dat     fort.96 
	    ln -s $bcphx/z_m0d5.afe_$bcsuf.dat     fort.97 
	    ln -s $bcphx/z_p0d0.afe_$bcsuf.dat     fort.98 
	    ln -s $bcphx/z_p0d3.afe_$bcsuf.dat     fort.99
	elif (( $fehint <= -75 )); then
	    echo "[Fe/H] = -1.0 BC used"
	    ln -s $bckur/atmk1990m10.tab           fort.38 
	    ln -s $bcphx/z_m1d0.afe_$bcsuf.dat     fort.95 
	    ln -s $bcphx/z_m0d5.afe_$bcsuf.dat     fort.96 
	    ln -s $bcphx/z_p0d0.afe_$bcsuf.dat     fort.97 
	    ln -s $bcphx/z_p0d3.afe_$bcsuf.dat     fort.98 
	    ln -s $bcphx/z_p0d5.afe_$bcsuf.dat     fort.99
	elif (( $fehint <= -25 )); then
	    echo "[Fe/H] = -0.5 BC used"
	    ln -s $bckur/atmk1990m05.tab           fort.38 
	    ln -s $bcphx/z_m1d0.afe_$bcsuf.dat     fort.95 
	    ln -s $bcphx/z_m0d5.afe_$bcsuf.dat     fort.96 
	    ln -s $bcphx/z_p0d0.afe_$bcsuf.dat     fort.97 
	    ln -s $bcphx/z_p0d3.afe_$bcsuf.dat     fort.98 
	    ln -s $bcphx/z_p0d5.afe_$bcsuf.dat     fort.99
	elif (( $fehint <= 15 )); then
	    echo "[Fe/H] = 0.0 BC used"
	    ln -s $bckur/atmk1990p00.tab           fort.38 
	    ln -s $bcphx/z_m1d0.afe_$bcsuf.dat     fort.95 
	    ln -s $bcphx/z_m0d5.afe_$bcsuf.dat     fort.96 
	    ln -s $bcphx/z_p0d0.afe_$bcsuf.dat     fort.97 
	    ln -s $bcphx/z_p0d3.afe_$bcsuf.dat     fort.98 
	    ln -s $bcphx/z_p0d5.afe_$bcsuf.dat     fort.99
	elif (( $fehint <= 40 )); then
	    echo "[Fe/H] = +0.3 BC used"
	    ln -s $bckur/atmk1990p03.tab           fort.38 
	    ln -s $bcphx/z_p0d3.afe_$bcsuf.dat     fort.95 
	    ln -s $bcphx/z_p0d0.afe_$bcsuf.dat     fort.96 
	    ln -s $bcphx/z_m0d5.afe_$bcsuf.dat     fort.97 
	    ln -s $bcphx/z_m1d0.afe_$bcsuf.dat     fort.98 
	    ln -s $bcphx/z_m1d5.afe_$bcsuf.dat     fort.99
	else
	    echo "[Fe/H] = +0.5 BC used"
	    ln -s $bckur/atmk1990p05.tab           fort.38 
	    ln -s $bcphx/z_m1d0.afe_$bcsuf.dat     fort.95 
	    ln -s $bcphx/z_m0d5.afe_$bcsuf.dat     fort.96 
	    ln -s $bcphx/z_p0d0.afe_$bcsuf.dat     fort.97 
	    ln -s $bcphx/z_p0d3.afe_$bcsuf.dat     fort.98 
	    ln -s $bcphx/z_p0d5.afe_$bcsuf.dat     fort.99
	fi
    fi
}
