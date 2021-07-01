#
## Shaft for torque
#
algebraic3d

# shaft core

solid ShaftOut = cylinder ( 0, 0, -45; 0, 0, 45; 1.25)
	and plane (0, 0, -16; 0, 0, -1)
	and plane (0, 0,  16; 0, 0, 1);

solid ShaftIn = cylinder ( 0, 0, -50; 0, 0, 50; 0.125)
	and plane (0, 0, -17; 0, 0, -1)
	and plane (0, 0,  17; 0, 0, 1);

solid Shaft = ShaftOut and not ShaftIn -maxh=0.125;


# magnetizing coil

solid MagCoilOut = cylinder ( 0, 0, -45; 0, 0, 45; 2.5)
	and plane (0, 0, -16; 0, 0, -1)
	and plane (0, 0,  16; 0, 0, 1);
	
solid MagCoilIn = cylinder ( 0, 0, -50; 0, 0, 50; 1.75)
	and plane (0, 0, -17; 0, 0, -1)
	and plane (0, 0,  17; 0, 0, 1);

solid MagCoil = MagCoilOut and not MagCoilIn -maxh=0.75;




# sensing coil

solid SensCoilOut = cylinder (0, 0, -7; 0, 0, 7; 3.75)
	and plane (0, 0, -3; 0, 0, -1)
	and plane (0, 0,  3; 0, 0, 1);

solid SensCoilIn = cylinder (0, 0, -8; 0, 0, 8; 3)
	and plane (0, 0, -4; 0, 0, -1)
	and plane (0, 0,  4; 0, 0, 1);

solid SensCoil = SensCoilOut and not SensCoilIn -maxh=0.75;


# Sphere

solid Range = sphere (0, 0, 0; 160) 
	 and not Shaft
	 and not MagCoil
	 and not SensCoil; # -maxh=30;  # 0.2


tlo Shaft -col=[0,0,1];

tlo MagCoil -col=[0,1,0];

tlo SensCoil -col=[1,0,0];

tlo Range -col=[0.5,0,0.5] -transparent;
