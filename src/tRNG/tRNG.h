/*	tRNG library
	true Random Number Generator

	(P) 2004 CAFxX (www.cafxx.cjb.net)
	This library is free for non-commercial use.			*/

// Includes
#include <windows.h>
#include <winbase.h>
#include <stdlib.h>
#include <stddef.h>
#include <process.h>
#include <memory.h>
#include <mmsystem.h>
#include <stdio.h>
#include <crtdbg.h>
#include <math.h>

// Error codes used in Status
#define TRNG_NO_ERROR                   0
#define TRNG_SHUTDOWN				   -1
#define TRNG_BAD_ALLOCATION             1
#define TRNG_CANT_START_GENERATOR       2
#define TRNG_MMSYSTEM_ERROR				3
#define TRNG_MMSYSTEM_NO_SOUNDCARD		4
#define TRNG_MMSYSTEM_WRONG_FORMAT		5
#define TRNG_NOT_ENOUGH_ENTROPY			6

// Settings flags
#define TRNG_ENABLE_QUALITY_CONTROL		1 // RECOMENDED
#define TRNG_ENABLE_UNBIASER			2 // NOT WORKING
#define TRNG_ENABLE_FHSPRNG				4 // NOT YET PORTED TO LIBRARY
#define TRNG_ENABLE_AUTOMATIC_MODE	    8 
#define TRNG_DISABLE_LOOP_COUNT		   16 // RECOMENDED

// Default values
#define DEFAULT_BUFFER_SIZE         40000 // bytes
#define DEFAULT_ITERATIONS			    8
#define DEFAULT_LOOPS_LIMIT			   50
#define DEFAULT_WAIT_INTERVAL		   20 // milliseconds
#define DEFAULT_PRNG_INSTANCES			8
#define DEFAULT_PRIME_NUMBER			3

class tRNG {
    public:
							tRNG(int BufferSize, int Iterations, int Settings);
							~tRNG();
        int					rand();
        
        int					Error();
		static const int	Version = 3;
		void				Dump();

    private:
        void				Generate();
        void				CreateGeneratorThread();
		friend void			ExecuteGenerator(void*);
		void				InitSource(); 
		void				CloseSource();
		bool				ControlQuality();
		void				Unbias();
		int					EstimateEntropy(void*);

		int					prand();
		int					Seed[DEFAULT_PRNG_INSTANCES];

        char*				DispenserBuffer;
        char*				GeneratorBuffer;
        int					BufferSize;
		int					Iterations;
        int					DispenserBufferPosition;
        bool				GeneratorBufferIsReady;
		int					GeneratorThreadID;
        bool				GeneratorIsRunning;
		bool				GeneratorReset;
		WAVEINCAPS			waveincaps;
		HWAVEIN				waveinputhandle;      
		WAVEFORMATEX		waveinputformat;
		WAVEHDR				waveinputheader;

		int					Status;
		int					Underruns;
		int					BuffersGenerated;
		int					BytesServed;
		bool				UseUnbiasing;
		bool				UseQualityControl;
		bool				UseFHSpRNG;
		bool				UseAutomaticMode;
		bool				DisableLoop;
};  

const int SampleRateCount = 9;	const int SampleRateList[] = { 192000, 96000, 48000, 44100, 32000, 22050, 16000, 11025, 8000 };
const int ChannelCount = 2;		const int ChannelList[] = { 2, 1 };
const int SampleDepthCount = 2;	const int SampleDepthList[] = { 16, 8 }; // also 32 and 24 bps should be possible, but they don't work. (?)
