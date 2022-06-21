/*
 COMPILE with:
    gcc -o L6_ofdm_time L6_ofdm_timing_sync.c Transmitter.c Receiver.c Channel.c Auxiliary.c -lm -Wall

 Usage example:
    ./L6_ofdm_time qpsk mmse 1 5000

*/

#include "OFDM.h"


void usage(char* progName) {
    printf("\nUsage: %s <num_taps> <symbol blocks> <pilot rate>\n\n", progName);
    printf("modulation: [ BPSK | bpsk | QPSK | qpsk | 16QAM | 16qam | 64QAM | 64qam | 256QAM | 256qam ]\n");
    printf("equalizer:  [ ZF | zf | MRC | mrc | MMSE | mmse ]\n");
    printf("            MRC is not appropriate for QAM.\n");
    printf("num_taps: length of the channel delay profile\n");
    printf("symbol blocks: number of symbol blocks to simulate\n");
}

void output_file_name(char const* filePath, char* fileName, char const* varname,
                      char const* modType, char const* eqType, char const* numTaps) {

    // modulation and equalization type strings' lengths
    int lenMod = strlen(modType);
    int lenEq = strlen(eqType);
    
    // ensure modulation and equalizer name are in upper case
    char modName[10], eqName[10];

    strcpy(modName, modType);
    strcpy(eqName, eqType);
    
    for (int i = 0; i < lenMod; i++) {
        modName[i] = toupper( modType[i] );
    }
    for (int i = 0; i < lenEq; i++) {
        eqName[i] = toupper( eqType[i] );
    }

    // assemble the file name
    strcpy(fileName, filePath);
    strcat(fileName, varname);
    strcat(fileName, "_");

    strcat(fileName, modName);
    
    strcat(fileName, "_");
    strcat(fileName, eqName);

    strcat(fileName, "_L");
    strcat(fileName, numTaps);
    strcat(fileName, ".txt");
}


int main(int argc, char** argv) {
    
    // Initialize the random generator. Must be done once per run of main
    srand((unsigned) time(NULL));
    
    // set if SNR-per-bit values for BER/SER evalation
    double EbN0[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40}; 

    // check required if all required parameters are provided
    if (argc != 5) {
        usage(argv[0]);
        return(EXIT_FAILURE);
    }
    
    // set number of bits per symbol
    int bitsPerSymbol = get_bits_per_symbol(argv[1]);

    // check if apppropriate modulation was specified
    if (bitsPerSymbol == -1) {
        usage(argv[0]);
        return(EXIT_FAILURE);
    }

    // set equalizer type from input paremeter (report error if MRC is used with QAM)
    equalizerType equalizer;
    if ( set_equalizer_type(bitsPerSymbol, argv[2], &equalizer) == EXIT_FAILURE) { 
        return EXIT_FAILURE;
    }

    // read the number of channel taps from parameter list
    int numTaps = atoi(argv[3]);
    
    // set number of blocks to simulate
    int numBlocks = atoi(argv[4]);
    int bitsPerBlock = SYMBOLS_PER_BLOCK * bitsPerSymbol; // bits per block (OFDM symbol)

    // Number of OFDM symbols in each slot (Preamble + Data)
    int slotLen = 2;

    // Parameters for timing synchronization
    int burstLen = 40;          // maximum value for timing offset
    int timingOffsetValue;      // random timing offset
    int estimatedTimingOffset;  // estimated timing offset

    // Set parameters for m-sequence
    int seqLen = 127;           // the length of m-sequences
    int syncSigStartIndex = 55; // first subcarrier to put m and gold sequence
    int Nid_2 = 2;              // Cell ID 2

    // Throughput loss due to GI insertion
    double cpLoss = (double) SYMBOLS_PER_BLOCK / (SYMBOLS_PER_BLOCK + CP_LEN);

    // Signal magnitude adjustment to match specified EbNo
    // NOTE: Noise is assumed to be of unit power per I/Q component
    double* snr = (double*) malloc(sizeof(EbN0));
    double* sqrtSNR = (double*) malloc(sizeof(EbN0));
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {
        snr[i] = pow( 10, EbN0[i] / 10); // SNR per bit
        
        // Adjust SNR per symbol to match Eb/No
        snr[i] *= bitsPerSymbol;    // SNR per symbol
        snr[i] *= cpLoss;           // GI (cyclic prefix) insertion
        sqrtSNR[i] = sqrt(snr[i]);  // amplitude scaling
    }

    // arrays to keep Tx and Rx bits in each iteration (i.e. one OFDM symbol)
    unsigned char* txBits = (unsigned char*) malloc( bitsPerBlock * sizeof(unsigned char));
    unsigned char* rxBits = (unsigned char*) malloc( bitsPerBlock * sizeof(unsigned char));

    // Tx symbols
    double* txSymI = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));
    double* txSymQ = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));
    
    // Tx OFDM symbol (modulated)
    double* txModI = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* txModQ = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    
    // Rx OFDM symbol (after channel)
    double* rxSymI = (double*) calloc( slotLen* (SYMBOLS_PER_BLOCK + CP_LEN), sizeof(double));
    double* rxSymQ = (double*) calloc( slotLen* (SYMBOLS_PER_BLOCK + CP_LEN), sizeof(double));
    
    // Rx OFDM symbol (after timing offset)
    double* rxSymOffsetI = (double*) calloc( slotLen * (SYMBOLS_PER_BLOCK + CP_LEN) + burstLen, sizeof(double));
    double* rxSymOffsetQ = (double*) calloc( slotLen * (SYMBOLS_PER_BLOCK + CP_LEN) + burstLen, sizeof(double));

    // Rx symbols
    double* rxEstI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* rxEstQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    
    // multipath channel
    double* hi = (double*) malloc(numTaps * sizeof(double));
    double* hq = (double*) malloc(numTaps * sizeof(double));
    
    double* alpha = (double*) malloc(N * numTaps * sizeof(double));
    double* phi = (double*) malloc(N * numTaps * sizeof(double));

    // Channel transfer function estimate
    double* HestI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* HestQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));

    // Bit and symbol error calculation
    double* ber = (double*) malloc(sizeof(EbN0));

    // Preamble 
    // Generating preamble for synchronization
    double* preambleI = (double*) calloc(SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* preambleQ = (double*) calloc(SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));

    // Tx preamble symbol
    double* txPreambleModI = (double*) calloc(SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* txPreambleModQ = (double*) calloc(SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));

    // Tx slot
    double* txSlotI = (double*) calloc(slotLen * (SYMBOLS_PER_BLOCK + CP_LEN), sizeof(double));
    double* txSlotQ = (double*) calloc(slotLen * (SYMBOLS_PER_BLOCK + CP_LEN), sizeof(double));

    // Preamble Symbol
    m_seq(preambleI + syncSigStartIndex + CP_LEN, preambleQ + syncSigStartIndex + CP_LEN, seqLen, Nid_2);
    ifft(preambleI + CP_LEN, preambleQ + CP_LEN, SYMBOLS_PER_BLOCK, txPreambleModI + CP_LEN, txPreambleModQ + CP_LEN);
    insert_cp( SYMBOLS_PER_BLOCK, CP_LEN, txPreambleModI, txPreambleModQ);
    
    // >> SIMULATION <<

    // terminal output header
    printf("\nEbNo\tBER\n");

    // Initialize/reset multipath fading simulator
    multipath_fading_init(N, numTaps, fDT, alpha, phi);

    // Repeat simulation for each given SNR point
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {
        
        // initialize bit/block error counters
        ber[i] = 0;

        // OFDM simulation loop (block by block)
        for (int j = 0; j < numBlocks; j++) {

            ///////////////////
            //  TRANSMITTER  //
            ///////////////////


            // generate data bits
            generate_info_bits(bitsPerBlock, txBits);
            
            // modulate bit sequence
            generate_symbols(txBits, bitsPerSymbol, SYMBOLS_PER_BLOCK, txSymI, txSymQ);

            // OFDM modulation
            ifft(txSymI, txSymQ, SYMBOLS_PER_BLOCK, 
                    txModI + CP_LEN, txModQ + CP_LEN);
            
            // insert cyclic prefix
            insert_cp(SYMBOLS_PER_BLOCK, CP_LEN, txModI, txModQ);
            
            //make a slot
            for (int i = 0; i < SYMBOLS_PER_BLOCK + CP_LEN; i++){
                txSlotI[i] = txPreambleModI[i];
                txSlotQ[i] = txPreambleModQ[i];
                txSlotI[i + SYMBOLS_PER_BLOCK + CP_LEN] = txModI[i];
                txSlotQ[i + SYMBOLS_PER_BLOCK + CP_LEN] = txModQ[i];
            }

            // scale Tx symbols to match SNR                    
            path_loss(txSlotI, txSlotQ, slotLen * (SYMBOLS_PER_BLOCK + CP_LEN), sqrtSNR[i]);
            
            ///////////////
            //  CHANNEL  //
            ///////////////

            // MULTIPATH FADING
            // update channel impulse response at the beginning of each block
            multipath_fading(alpha, phi, N, numTaps, j * (SYMBOLS_PER_BLOCK + CP_LEN), hi, hq);  // 
            
            // Convolution with channel impulse response
            conv_with_gi(hi, hq, numTaps, txSlotI, txSlotQ, slotLen * (SYMBOLS_PER_BLOCK + CP_LEN), 
                    0, rxSymI, rxSymQ);
            
            // TIMING OFFSET
            // generate random timing offset for each block and apply it to each block
            timingOffsetValue = rand() % burstLen;

            for (int i = 0; i < slotLen * (SYMBOLS_PER_BLOCK + CP_LEN); i++) {
                rxSymOffsetI[i + timingOffsetValue] = rxSymI[i];
                rxSymOffsetQ[i + timingOffsetValue] = rxSymQ[i];
            }

            // AWGN
            add_noise(rxSymOffsetI, rxSymOffsetQ, slotLen * (SYMBOLS_PER_BLOCK + CP_LEN) + burstLen);
            
            ////////////////
            //  RECEIVER  //
            ////////////////

            // TIMING SYNCHRONIZATION
            estimatedTimingOffset = timing_synch(rxSymOffsetI, rxSymOffsetQ, slotLen * (SYMBOLS_PER_BLOCK + CP_LEN) + burstLen);

            // Automatic Gain Control (AGC)
            // SNR estimation is typically done based on preamble (perfect estimate is assumed here)
            // NOTE: Reverts the pilot scaling introduced in path_loss function
            for (int k = estimatedTimingOffset + 2 * CP_LEN + SYMBOLS_PER_BLOCK; k < estimatedTimingOffset + SYMBOLS_PER_BLOCK + 2 * CP_LEN + SYMBOLS_PER_BLOCK; k++) {
                rxSymOffsetI[k] /= (sqrtSNR[i] * SQRT_OF_2);
                rxSymOffsetQ[k] /= (sqrtSNR[i] * SQRT_OF_2);
            }

            // OFDM demodulation (CI is dropped)
            fft(rxSymOffsetI + 2 * CP_LEN + SYMBOLS_PER_BLOCK + estimatedTimingOffset, rxSymOffsetQ + 2 * CP_LEN + SYMBOLS_PER_BLOCK + estimatedTimingOffset, SYMBOLS_PER_BLOCK, rxEstI, rxEstQ);
            

            // Perfect channel knowledge
            for (int k = 0; k < SYMBOLS_PER_BLOCK; k++) {
                HestI[k] = 0;
                HestQ[k] = 0;
                for (int m = 0; m < numTaps; m++) {
                    HestI[k] += hi[m] * cos(-2 * PI * k * m / SYMBOLS_PER_BLOCK) 
                            - hq[m] * sin(-2 * PI * k * m / SYMBOLS_PER_BLOCK);
                    HestQ[k] += hq[m] * cos(-2 * PI * k * m / SYMBOLS_PER_BLOCK) 
                            + hi[m] * sin(-2 * PI * k * m / SYMBOLS_PER_BLOCK);
                }
            }

            // Equalization (takes channel transfer function as input)
            fde(equalizer, snr[i], HestI, HestQ, rxEstI, rxEstQ, SYMBOLS_PER_BLOCK);
            
            // Demodulate (detect symbols and output data bits)
            decode_symbols(rxEstI, rxEstQ, SYMBOLS_PER_BLOCK, bitsPerSymbol, rxBits);
            
            // Count bit and symbol errors in the block (only for data symbols)
            for (int m = 0; m < bitsPerBlock; m++) {
                if (txBits[m] != rxBits[m])  ber[i]++;
            }
        }
        
        // Calculate BER
        ber[i] /= (bitsPerBlock * numBlocks);
        printf("%f\t%f\n", EbN0[i], ber[i]);
    }
    
    // generate output file name
    char fileNameBER[FILE_NAME_SIZE];
    output_file_name("./results/", fileNameBER, "BER_TimingOffset", argv[1], argv[2], argv[3]);
    
    // save BER to file
    save_asignal(EbN0, ber, sizeof(EbN0)/sizeof(double), fileNameBER);
    printf("\nBER saved to %s\n", fileNameBER);
    
    // Release all allocated memory
    free(txBits);
    free(rxBits);
    free(preambleI);
    free(preambleQ);
    free(txSymI);
    free(txSymQ);
    free(txModI);
    free(txModQ);
    free(rxSymI);
    free(rxSymQ);
    free(rxSymOffsetI);
    free(rxSymOffsetQ);
    free(txPreambleModI);
    free(txPreambleModQ);
    free(txSlotI);
    free(txSlotQ);
    free(rxEstI);
    free(rxEstQ);
    free(snr);
    free(sqrtSNR);
    free(alpha);
    free(phi);
    free(hi);
    free(hq);
    free(ber);
    //
    free(HestI);
    free(HestQ);
    //
    return EXIT_SUCCESS;
}

