#include "OFDM.h"

/* Get power allocation across subcarriers according to the
 * WATERFILLING algorithm.
 *
 * Params
 *  powAlloc    - power allocation array (output)
 *  gamma       - channel SNR, i.e. |H(k)|^2/No
 *  numCarriers - numumber of subarriers, i.e. length of gamma and powAlloc 
 *  sqrtSNR     - square root of SNR
 * 
 * Output
 *  int - number of active subcarriers 
 * 
 */
int power_allocation(double powAlloc[], double gamma[], int numCarriers, double sqrtSNR) {

    //***Student_code_start*** 
	double Cth = 0.0, sum_ck = 0.0, count = 0;
	for(int m = 0; m < numCarriers; m++){
		sum_ck += 1/gamma[m]; 
	}
	Cth = numCarriers/(sqrtSNR*sqrtSNR*numCarriers*2+sum_ck);
	for(int i = 0; i < numCarriers; i++){
		powAlloc[i] = 0;
		for(int j = 0; j < numCarriers; j++){
			powAlloc[i] = (1/Cth-1/gamma[j]) > 0?(1/Cth-1/gamma[j]):0;
			
		}
		if(powAlloc[i] != 0){count++;}
	}
	return count;
	//***Student_code_end*****

}


/* Allocate bits according to the SNR, by maximizing the constellation order
 * while ensuring that the target BER is met.
 *
 * Params
 *  bitAlloc       - bit allocation array (output)
 *  gamma          - channel SNR, i.e. |H(k)|^2/No > assumes E{|H(k)|^2} = 1
 *  numCarriers    - number of subcarriers, i.e. the length of bitAlloc, powAlloc and gamma 
 *  snrAvg         - average Rx SNR (estimated at the Rx)
 *  snrTh          - SNR thresholds for modulation up to the largest order
 *  modsBitsPerSym - array with number of bits per symbol for available modulations
 *  numMods        - number of considered modulations, i.e. length of modsBitsPerSym
 * 
 * Output
 *  int - total number of allocated bits
 * 
 */
int bit_allocation(int bitAlloc[], double gamma[], int numCarriers,
                 double snrAvg, int modBitsPerSym[], double snrTh[], int numMods) {

    //***Student_code_start*** 
	int n, total_bits = 0;
	
	for(int i = 0; i < numCarriers; i++){
		n = 0;
		for(int j = 1; j < numMods; j++){
			if(snrAvg*gamma[i] > snrTh[j]){n+=1;}
		}
		//printf("%d",i);
		//printf("%d",n);
		bitAlloc[i] = modBitsPerSym[n];
		total_bits += bitAlloc[i];
	}
	//***Student_code_end*****
	return total_bits;
}


/* Modulate sequence with variable modulation order (rate adaptive case).
 *
 * Parameters
 *  txBits          -
 *  bitsPerSymbol   -
 *  numSymbols      -
 *  txSymI          -
 *  txSymQ          -
 */
void generate_symbols_ra(unsigned char txBits[], int bitsPerSymbol[], int numSymbols,
        double txSymI[], double txSymQ[]) {
    
    //***Student_code_start*** 
	int count = 0;
	for(int i = 0; i < numSymbols; i++){
		generate_asymbol(txBits+count, bitsPerSymbol[i], txSymI+i, txSymQ+i);
		count += bitsPerSymbol[i];
	}
	//***Student_code_end*****

}

/* Demodulate sequence with variable modulation order (rate adaptive case).
 *
 * Parameters
 *  rxSymI          -
 *  rxSymQ          -
 *  numSymbols      -
 *  bitsPerSymbol   -
 *  rxBits          -
 */
void decode_symbols_ra(double rxSymI[], double rxSymQ[], int numSymbols,
        int bitsPerSymbol[], unsigned char rxBits[]) {
    
    //***Student_code_start*** 
	int count = 0;
	for(int i = 0; i < numSymbols; i++){
		decode_asymbol(rxSymI+i, rxSymQ+i, rxBits+count, bitsPerSymbol[i]);
		count += bitsPerSymbol[i];
	}
	//***Student_code_end*****

}