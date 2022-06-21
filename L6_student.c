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

    //***Student_code_start*** 
	//    Your code...
	//    ...
	//    ...
	//***Student_code_end*****

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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
	//***Student_code_end*****
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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
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
	//    Your code...
	//    ...
	//    ...
	//***Student_code_end*****
    
}


