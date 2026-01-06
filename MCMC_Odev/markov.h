#ifndef __MARKOV_H__
#define __MARKOV_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <time.h>




class MARKOV
{
	public:
		MARKOV();
		void estimateNucleotidsInverseTransform();


	private:
		void 	reinitializeParams();
		double 	get_random(void);
		int 	it_sampler(int NACode);
		void 	calculate_transition_matrix(const char* array,uint32_t length);


		//Nucleic Acids
		const int 	m_NumberOfNA;
		char 		m_NucleicAcids[4];

		double 		m_TrMatrix[4][4]; 	
		double 		m_NAHist[4];	  	
		double 		m_StateHist[4][4];
		// stationary probability distribution after monte carlo
		int 		StPrDist[4];
};












#endif
