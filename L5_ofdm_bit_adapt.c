/*
RATE ADAPTATION BER simulation

 COMPILE with:
    gcc -o L5_ofdm_bit_adapt L5_ofdm_bit_adapt.c Transmitter.c Receiver.c Channel.c Auxiliary.c -lm -Wall

 Usage example:
    ./L5_ofdm_bit_adapt 8 5000

 Zip folder on server for download (call from 'Documents' folder):
    zip -r L5_new.zip L5_new
*/

#include "OFDM.h"


void usage(char* progName) {
    printf("\nUsage: %s <num_taps> <symbol blocks> <pilot rate>\n\n", progName);
    printf("num_taps: length of the channel delay profile\n");
    printf("symbol blocks: number of symbol blocks to simulate\n");
}


void output_file_name(char const* filePath, char const* varname, char* fileName, char* howManyTaps) {

    strcpy(fileName, filePath);
    strcat(fileName, varname);
    strcat(fileName, "_RA_MMSE"); // MMSE equalizer
    // strcat(fileName, "_WF_ZF"); // ZF equalizer

    strcat(fileName, "_L");
    strcat(fileName, howManyTaps);
    strcat(fileName, ".txt");
}

int main(int argc, char** argv) {
    
    // Initialize the random generator. Must be done once per run of main
    srand((unsigned) time(NULL));
    
    // set if SNR-per-bit values for BER/SER evalation
    double snr_dB[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50};


    // check required if all required parameters are provided
    if (argc != 3) {
        usage(argv[0]);
        return(EXIT_FAILURE);
    }

    // set equalizer type
    equalizerType equalizer = MMSE;
    // equalizerType equalizer = ZF;
    
    // read the number of channel taps from parameter list
    int numTaps = atoi(argv[1]);
    
    // set number of blocks to simulate
    int numBlocks = atoi(argv[2]);


    // considered modulations (i.e. corresponding bits per symbol)
    int modBitsPerSym[] = {0, QPSK, QAM16, QAM64, QAM256};

    // number of consider modulations
    int numMods = sizeof(modBitsPerSym)/sizeof(modBitsPerSym[0]);

    // target bit error rate
    double ber_target =  1e-3;

    // maximum number of bits per symbol (the largest constellation)
    int maxBitsPerSymbol = modBitsPerSym[numMods-1];

    // maximum number of bits per block (OFDM symbol)
    int maxBitsPerBlock = SYMBOLS_PER_BLOCK * maxBitsPerSymbol;

    // bit allocation array (num. bits per symbol for each subcarrier)
    int bitAlloc[SYMBOLS_PER_BLOCK];

    // num. of active subcarriers (all)
    int numActiveCarriers = SYMBOLS_PER_BLOCK;

    // number of bits per block (known after adaptive algorithm)
    int bitsPerBlock;
    
    // average number of bits per symbol
    double avgBitsPerSym;

    // arrays to keep SNR in linear domain (and its square root)
    double snr[sizeof(snr_dB)/sizeof(double)];
    double sqrtSNR[sizeof(snr_dB)/sizeof(double)];


    // arrays to keep Tx and Rx bits in each iteration (i.e. one OFDM symbol)
    unsigned char* txBits = (unsigned char*) malloc( maxBitsPerBlock * sizeof(unsigned char));
    unsigned char* rxBits = (unsigned char*) malloc( maxBitsPerBlock * sizeof(unsigned char));

    // Tx symbols
    double* txSymI = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));
    double* txSymQ = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));
    
    // Tx OFDM symbol (modulated)
    double* txModI = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* txModQ = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    
    // Rx OFDM symbol (after channel)
    double* rxSymI = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* rxSymQ = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    
    // Rx symbols
    double* rxEstI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* rxEstQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    
    // multipath channel
    double* hi = (double*) malloc(numTaps * sizeof(double));
    double* hq = (double*) malloc(numTaps * sizeof(double));
    
    double* alpha = (double*) malloc(NUM_RAYS * numTaps * sizeof(double));
    double* phi = (double*) malloc(NUM_RAYS * numTaps * sizeof(double));

    // Channel transfer function estimate
    double* HestI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* HestQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));

    
    // bit error calculation
    double* ber = (double*) malloc(sizeof(snr_dB));

    // average throughput calculation
    double* avgRate = (double*) malloc(sizeof(snr_dB));

    // total number of simulated bits (for a given EbNo point)
    int numBitsTot;
 
    // total number of simulated symbols (for a given EbNo point)
    int numSymbolsTot;

    // number of blocks with all subcarriers inactive
    int numEmptyBlocks;

    // array to keep channel SNR, i.e. gamma(f) = |H(f)|^2/No)
    double* gamma = (double*) malloc(SYMBOLS_PER_BLOCK * sizeof(double));

    // get the SNR thresholds for all modulation sizes (given the taget BER)
    double* snr_th = (double*) malloc( numMods * sizeof(double) );
    for(int k = 0; k < numMods; k++) {
        snr_th[k] = get_snr_threshold_alt(modBitsPerSym[k], ber_target); // Goldsmith BER approximation
    }
    

    // >> SIMULATION <<
    printf("\n\n");
    printf("LEGEND\n");
    printf("  SNR\t  Signal to Noise ratio\n");
    printf("  BER\t  Bit Error Probability\n");
    printf("  R_avg\t  Average throughput (in bits per block)\n");
    printf("  n_avg\t  Average number of bits per symbol\n");
    printf("  #Empty  Number of empty blocks (with no allocated bits)\n");
    printf("  n_tot   Total number of simulated bits\n");
    printf("\n");


    // terminal output header
    printf("\nSNR\tBER\t\tR_avg\tn_avg\t#Empty\tn_tot\n");

    // Repeat simulation for each given SNR point
    for (int i = 0; i < sizeof(snr_dB) / sizeof(double); i++) {

        // SNR 
        snr[i] = pow( 10, snr_dB[i] / 10);
        sqrtSNR[i] = sqrt(snr[i]);  // amplitude scaling

        
        // initialize counters
        numEmptyBlocks = 0;
        numSymbolsTot = 0;
        numBitsTot = 0;
        ber[i] = 0;
        avgRate[i] = 0;

        // Initialize/reset multipath fading simulator
        multipath_fading_init(NUM_RAYS, numTaps, fDT, alpha, phi);


        // OFDM simulation loop (block by block)
        for (int j = 0; j < numBlocks; j++) {


            ////////////////////////////////////////
            //  CHANNEL KNOWLEDGE AT TRANSMITTER  //
            ////////////////////////////////////////
            /* NOTE:
               Channel is now updated before the Tx processing, as the channel transfer
               function is required bit allocation algorithm.
               This emulates channel knowledge at the Tx, which is in practice obtained
               through feedback from the Rx. 
            */
            
            // uncorrelated fading
            multipath_fading_init(NUM_RAYS, numTaps, fDT, alpha, phi);
            multipath_fading(alpha, phi, NUM_RAYS, numTaps, 0, hi, hq);
            
            // Channel Transfer Function (Perfect channel knowledge)
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


            ///////////////////
            //  TRANSMITTER  //
            ///////////////////

            // BIT ALLOCATION

            // channel power transmission coefficient, |H(k)|^2
            for (int k = 0; k < SYMBOLS_PER_BLOCK; k++) {
                gamma[k] = HestI[k]*HestI[k] + HestQ[k]*HestQ[k];
            }

            // bit allocation (adaptive modulation)
            bitsPerBlock = bit_allocation(bitAlloc, gamma, SYMBOLS_PER_BLOCK, snr[i], modBitsPerSym, snr_th, numMods);

            // >>
            // if no bits are allocated, skip to the next frame
            if (bitsPerBlock == 0) {
                numEmptyBlocks++;
                continue;
            }

            // generate data bits
            generate_info_bits(bitsPerBlock, txBits);
            
            // modulate bit sequence with variable modulation order (bit loading)
            generate_symbols_ra(txBits, bitAlloc, SYMBOLS_PER_BLOCK, txSymI, txSymQ);


            // OFDM modulation
            ifft(txSymI, txSymQ, SYMBOLS_PER_BLOCK, 
                    txModI + CP_LEN, txModQ + CP_LEN);
            
            // insert cyclic prefix
            insert_cp(SYMBOLS_PER_BLOCK, CP_LEN, txModI, txModQ);

            // scale Tx symbols to match SNR                    
            path_loss(txModI, txModQ, SYMBOLS_PER_BLOCK + CP_LEN, sqrtSNR[i]);
            

            ///////////////
            //  CHANNEL  //
            ///////////////
            
            // Convolution with channel impulse response
            conv_with_gi(hi, hq, numTaps, txModI, txModQ, SYMBOLS_PER_BLOCK + CP_LEN, 
                    CP_LEN, rxSymI, rxSymQ);
            
            // AWGN
            // NOTE: Noise is assumed to be of unit power per I/Q component (total power is 2)
            add_noise(rxSymI, rxSymQ, SYMBOLS_PER_BLOCK + CP_LEN);
            

            ////////////////
            //  RECEIVER  //
            ////////////////

            // Automatic Gain Control (AGC)
            // SNR estimation is typically done based on preamble (perfect estimate is assumed here)
            // NOTE: Reverts the pilot scaling introduced in path_loss function
            for (int k = CP_LEN; k < SYMBOLS_PER_BLOCK+CP_LEN; k++) {
                rxSymI[k] /= (sqrtSNR[i]*SQRT_OF_2);
                rxSymQ[k] /= (sqrtSNR[i]*SQRT_OF_2);
            }

            // OFDM demodulation (CI is dropped)
            fft(rxSymI + CP_LEN, rxSymQ + CP_LEN, SYMBOLS_PER_BLOCK, rxEstI, rxEstQ);
            
            // Equalization (takes channel transfer function as input)
            fde(equalizer, snr[i], HestI, HestQ, rxEstI, rxEstQ, SYMBOLS_PER_BLOCK);
            
            // Demodulate with variable modulation orders
            decode_symbols_ra(rxEstI, rxEstQ, SYMBOLS_PER_BLOCK, bitAlloc, rxBits);
            
            // Count bit errors in the block
            for (int k = 0; k < bitsPerBlock; k++) {
                if (rxBits[k] != txBits[k])  ber[i]++;
            }

            // increment number of transmitted (info) symbols and  bits
            numSymbolsTot += numActiveCarriers; // number of actual symbols in block
            numBitsTot += bitsPerBlock;         // for BER and throughput calculation
        }

        // average number of bits per symbol for the current SNR point
        avgBitsPerSym = ((double) numBitsTot)/((double)numSymbolsTot);

        // Calculate BER and average throughput
        ber[i] /= ((double) numBitsTot);
        avgRate[i] = ((double) numBitsTot)/((double) numBlocks - numEmptyBlocks);

        // print the current simulation loop state
        printf("%.0f\t%f\t%.2f\t%.2f\t%d\t%d\n", snr_dB[i], ber[i], avgRate[i], avgBitsPerSym, numEmptyBlocks, numBitsTot);
    }

    // generate output file name for BER
    char fileNameBER[FILE_NAME_SIZE];
    output_file_name("./results/", "BER", fileNameBER, argv[1]);
    
    // save BER to file
    save_asignal(snr_dB, ber, sizeof(snr_dB)/sizeof(double), fileNameBER);
    printf("\nBER saved to %s\n", fileNameBER);


    // generate output file name for throughput
    char fileNameRATE[FILE_NAME_SIZE];
    output_file_name("./results/", "RATE", fileNameRATE, argv[1]);

    // save throughput to file
    save_asignal(snr_dB, avgRate, sizeof(snr_dB)/sizeof(double), fileNameRATE);
    printf("RATE saved to %s\n\n", fileNameRATE);
    
    // Release all allocated memory
    free(txBits);
    free(rxBits);
    free(txSymI);
    free(txSymQ);
    free(txModI);
    free(txModQ);
    free(rxSymI);
    free(rxSymQ);
    free(rxEstI);
    free(rxEstQ);
    free(alpha);
    free(phi);
    free(hi);
    free(hq);
    free(HestI);
    free(HestQ);
    //
    free(ber);
    free(avgRate);
    //
    free(gamma);
    free(snr_th);
    //
    return EXIT_SUCCESS;
}

