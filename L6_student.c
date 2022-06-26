#include "OFDM.h"

/* m-sequence is written based on the formula
*
* preambleI   - I-component part of generated m-sequence
* preambleQ   - Q-component part of generated m-sequence
* len         - Length of the m-sequence
* Nid_2       - Value of second Cell ID 
*
*/ 
void m_seq(double preambleI[], double preambleQ[], int len, int Nid_2) {
	unsigned char* x = (unsigned char*) malloc((len) * sizeof(unsigned char));
	int m = 0;

	x[0] = 0;
	x[1] = 1;
	x[2] = 1;
	x[3] = 0;
	x[4] = 1;
	x[5] = 1;
	x[6] = 1;
	
	for(int i = 0; i < len-7; i++){
		x[i+7] = (x[i+4]+x[i])%2;
	}

	for(int n = 0; n < len; n++){
		m = (n+43*Nid_2) % 127;
		preambleI[n] = 1.0-2.0*x[m];
		preambleQ[n] = 0.0;
	}

}


/* Estimating timing offset value by autocorrelation on CP.
* 
* inI   - I-component of the received signal
* inQ   - Q-component of the received signal
* lenIN - Length of received signal
* outI  - I-component of the time synchronized signal
* outQ  - Q-component of the time synchronized signal
* normalizedValue  - Normalized function for time offset estimation
*
*/ 
void auto_corr(double* inI, double* inQ, int lenIN, double* outI, double* outQ, double* normalizedValue) {  
    
    //***Student_code_start*** 
	for(int l = 0; l < lenIN-2*SYMBOLS_PER_BLOCK-2*CP_LEN; l++){
		for(int i = 0; i < CP_LEN; i++){
			normalizedValue[l] += inI[l+i+SYMBOLS_PER_BLOCK]*inI[l+i+SYMBOLS_PER_BLOCK]+inQ[l+i+SYMBOLS_PER_BLOCK]*inQ[l+i+SYMBOLS_PER_BLOCK];
			outI[l] += inI[l+i]*inI[l+i+SYMBOLS_PER_BLOCK]+inQ[l+i]*inQ[l+i+SYMBOLS_PER_BLOCK];
			outQ[l] += inQ[l+i]*inI[l+i+SYMBOLS_PER_BLOCK]-inI[l+i]*inQ[l+i+SYMBOLS_PER_BLOCK];
		}
		normalizedValue[l] = (outI[l]*outI[l]+outQ[l]*outQ[l])/normalizedValue[l];
	}
	//***Student_code_end*****
}


/* Finding maximum value of an array
*
* inputArray    - Input array
* len           - length of the input array
*
* Output
* index         - index of the maximum value
*/
int find_max_index(double* inputArray, int len) {

    //***Student_code_start*** 
	int max_idx = 0;
	for(int i = 0; i < len; i++){
		if(inputArray[i]>inputArray[max_idx]){
			max_idx = i;
		}
	}
	return max_idx;
	//***Student_code_end*****

}

/* Timing synchronization
*
* inI    - I-component of the received signal
* inQ    - Q-component of the received signal
* lenIN  - Length of received signal
* 
* Output
* maxIndex  - Index of maximum value
*
*/
int timing_synch(double* inI, double* inQ, int lenIN){

    //***Student_code_start***
	int offset = 0;  
	double* outI = (double*) calloc(lenIN-2*SYMBOLS_PER_BLOCK-2*CP_LEN, sizeof(double));
	double* outQ = (double*) calloc(lenIN-2*SYMBOLS_PER_BLOCK-2*CP_LEN, sizeof(double));
	double* normalizedValue = (double*) calloc(lenIN-2*SYMBOLS_PER_BLOCK-2*CP_LEN, sizeof(double));

	auto_corr(inI, inQ, lenIN, outI, outQ, normalizedValue);
	offset = find_max_index(normalizedValue, lenIN-2*SYMBOLS_PER_BLOCK-2*CP_LEN);

	free(outI);
	free(outQ);
	free(normalizedValue);

	return offset;
	//***Student_code_end*****

}

/* Esimating fractional frequency offset
*
* inI     - I-componnent of the received signal 
* inQ     - Q-componnent of the received signal 
* lenIN            - Length of the input
*
* Output
* fracfreqOffset   - fractional frequency offset value
*/

double carrier_freq_estimation_fractional(double* inI, double* inQ, int lenIN){
   
    //***Student_code_start*** 
	double c = 0;
	double freq_offset_frac = 0;
	double outI = 0;
	double outQ = 0;

	for(int i = 0; i < CP_LEN; i++){
		outI += inI[i]*inI[i+SYMBOLS_PER_BLOCK]+inQ[i]*inQ[i+SYMBOLS_PER_BLOCK];
		outQ += inQ[i]*inI[i+SYMBOLS_PER_BLOCK]-inI[i]*inQ[i+SYMBOLS_PER_BLOCK];
	}

	freq_offset_frac = -0.5/PI*(atan(outQ/outI+c));

	return freq_offset_frac;
	//***Student_code_end***
}

/* Finding index of maximum value of complex array
* 
* inI               - I-component of the received signal
* inQ               - Q-component of the received signal
* lenIN             - Length of received signal
*
* Output
* convIndex  - Index of maximum value 
*
*/
int find_peak_index(double* inI, double* inQ, int lenIN) {
    
    //***Student_code_start*** 
	int peak_idx = 0;
	double max = 0;

	for(int i = 0; i < lenIN; i++){
		if(inI[i]*inI[i]+inQ[i]*inQ[i] > max){
			max = inI[i]*inI[i]+inQ[i]*inQ[i];
			peak_idx = i;
		}
	}

	return peak_idx;
	//***Student_code_end*****
}

/* Estimating integer frequency offset
*
* rxSymI     - I-component of the received signal
* rxSymQ     - Q-component of the received signal
* seqLen     - length of m sequences
* syncSigStartIndex   - First subcarrier to put m-sequence
* Nid_2      - Cell ID 2
*
* Output
* carrierFreqOffsetInt   - Estimated value for integer carrier frequency offset
*/
int carrier_freq_estimation_integer(double* rxSymI, double* rxSymQ, int seqLen, int syncSigStartIndex, int Nid_2){
    
    //***Student_code_start*** 
	int freq_offset_int = 0;
	double* pI = (double*) calloc(seqLen, sizeof(double));
	double* pQ = (double*) calloc(seqLen, sizeof(double));
	double* rI = (double*) calloc(SYMBOLS_PER_BLOCK, sizeof(double));
	double* rQ = (double*) calloc(SYMBOLS_PER_BLOCK, sizeof(double));
	double* outI = (double*) calloc(SYMBOLS_PER_BLOCK-seqLen+1, sizeof(double));
	double* outQ = (double*) calloc(SYMBOLS_PER_BLOCK-seqLen+1, sizeof(double));
	
	m_seq(pI, pQ, seqLen, Nid_2);
	fft(rxSymI+CP_LEN, rxSymQ+CP_LEN, SYMBOLS_PER_BLOCK, rI, rQ);

	cross_corr(pI, pQ, seqLen, rI, rQ, SYMBOLS_PER_BLOCK, outI, outQ, SYMBOLS_PER_BLOCK-seqLen+1);
	freq_offset_int = find_peak_index(outI, outQ, SYMBOLS_PER_BLOCK-seqLen);
	freq_offset_int -= syncSigStartIndex;

	free(pI);
	free(pQ);
	free(rI);
	free(rQ);
	free(outI);
	free(outQ);
	
	return freq_offset_int;
	//***Student_code_end*****

}

/* Cross-Correlation of two complexsignals
*
* hI     - I-component of the first signal
* hQ     - Q-component of the first signal
* lenH   - length of the first signal
* inI    - I-component of the second signal
* inQ    - Q-component of the second signal
* lenIN  - length of the second signal
* outI   - I-component of the correlation output
* outQ   - Q-component of the correlation output
* lenOUT - length of the correlation output
*
*/
void cross_corr(double* hI, double* hQ, int lenH, double* inI, double* inQ, 
        int lenIN, double* outI, double* outQ, int lenOUT) {

    //***Student_code_start*** 
	for(int l = 0; l < lenOUT; l++){
		for(int i = 0; i < lenH; i++){
			outI[l] += inI[i+l]*hI[i] + inQ[i+l]*hQ[i];
			outQ[l] += inQ[i+l]*hI[i] - inI[i+l]*hQ[i];
		}
	}
	//***Student_code_end*****
    
}


