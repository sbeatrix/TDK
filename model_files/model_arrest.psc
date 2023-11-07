# Reactions
R1:
	$pool > c
   	ksyn_c20

R2:
   	$pool > cycb
    	ksyn_cy
	
R3:
    	c > $pool
    	kdegbg*c
	
R4:
     	mstar > m
     	kinact*mstar/(J+mstar)
	
	
R5:
    	m  > mstar
        kact*cycb*m*k_nuk*nuk/((J1+nuk)*(J+m))

R6:
    	mstar + c > mc
    	kassmc*mstar*c
	
R7:
    	mc > mstar + c
    	kdissmc*mc

R8:
    	a + c > ac
    	kassac*a*c

R9:
    	mc + ac > acmc
    	kassacmc*ac*mc

R10:
    	ac > a + c
    	kdissac*ac

R11:
    	acmc > mc + ac
    	kdissacmc*acmc

R12:
    	acmc > a + c + m
    	kdeg*acmc

R13:
    	cycb + ac > ac
    	kdeg_cy*cycb*ac
	
R14:
	mc > m
	kdegbg*mc
	
R15:
	ac > a
	kdegbg*ac 

R16:
	acmc > a + m
	kdegbg*acmc
	
R17:
	cycb > $pool
	kdegbg_cy*cycb
R18:
	nuk > $pool
	nuk*knUK_att*P


# InitPar
ksyn_c20 = 36.97731502237753
ksyn_cy = 316.94570558508684
kdegbg = 0.022521686855703566
kinact = 0.45009544921120104
J = 0.06168169109783722
kact = 0.0014200271790012738
k_nuk = 78
J1 = 0.5
kassmc = 0.021973634967776146
kdissmc = 0.6256358653434726
kassac = 0.023088755284985777
kassacmc = 0.04477093007618703
kdissac = 0.664021730464077
kdissacmc = 0.43441487879976926
kdeg = 1.117610105994406
kdeg_cy = 0.023952471280495782
kdegbg_cy = 0.0677970961252346
nuk=10
knUK_att=1
P=0

# InitVar
c = 25
mstar = 106
mc = 38
ac = 30
acmc = 30
m = 1
a = 80
cycb = 430



# Assignment rules
!F ctot = c + mc + ac + 2*acmc

FIX: P
Event: start_attachment, _TIME_ >=18000.0, 0.0


{
P=1
}
