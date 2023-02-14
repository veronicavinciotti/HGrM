#include <windows.h>
#include <winbase.h>
#include <stdlib.h>
#include <stddef.h>
#include <process.h>
#include <memory.h>
#include <mmsystem.h>

#define TRNG_NO_ERROR                   0
#define TRNG_BAD_ALLOCATION             1
#define TRNG_CANT_START_GENERATOR       2
#define TRNG_MMSYSTEM_ERROR				3
#define TRNG_MMSYSTEM_NO_SOUNDCARD		4

/*
class tRNG {
    public:
              tRNG(int BufferSize, int Iterations);
              ~tRNG();
        int   rand();
        
        int   Status;
    private:
        void  Generate(); // generator function to be executed asinchronously
        void  CreateGeneratorThread(); // function that executes Generate
		friend void ExecuteGenerator(void*);
        
        char* DispenserBuffer;
        char* GeneratorBuffer;
        int   BufferSize;
		int   Iterations;
        int   DispenserBufferPosition;
        bool  GeneratorBufferIsReady;
		int   GeneratorThreadID;
        bool  GeneratorIsRunning;

		// MMSYSTEM functions
		void  InitSource(); // initialize the soundcard
		void  RecordSource(); // start recording
		void  StopSource(); // stop recording
		void  CloseSource(); // close the soundcard

		WAVEINCAPS waveincaps;
		HWAVEIN waveinputhandle;      
		WAVEFORMATEX waveinputformat;
		WAVEHDR waveinputheader;
};  

*/