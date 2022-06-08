Bit:
		gcc -std=c99 -o L5_ofdm_bit_adapt L5_ofdm_bit_adapt.c L5_student.c -I include lib/libOFDM.so -lm -Wall

Power: 
	 gcc -std=c99 -o L5_ofdm_power_adapt L5_ofdm_power_adapt.c L5_student.c -I include lib/libOFDM.so -lm -Wall
	


