Freq:
		gcc -std=c99 -o L6_ofdm_freq_sync L6_ofdm_freq_sync.c L6_student.c -I include lib/libOFDM.so -lm -Wall

time: 
	 gcc -std=c99 -o L6_ofdm_timing_sync L6_ofdm_timing_sync.c L6_student.c -I include lib/libOFDM.so -lm -Wall
	


