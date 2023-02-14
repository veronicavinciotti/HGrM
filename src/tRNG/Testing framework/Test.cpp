#include "..\tRNG.h"
#include <stdio.h>
#include <time.h>

int main() {
	FILE* out = fopen("out.dat", "wb");
	int count = 40000;
	int* array = (int*)malloc(count*4);
	clock_t start, start2, end, end2;

	printf("tRNG library test framework\n");
	start = clock();
	{
    tRNG random(0,0,TRNG_ENABLE_QUALITY_CONTROL | TRNG_DISABLE_LOOP_COUNT);
	start2 = clock();
    for (int i = 0; i < count; i++) {
		array[i] = random.rand();
//        printf("%08X", random.rand()); 
    }
	end = clock();
	fwrite(array, i, 4, out);
	}
	end2 = clock();
	printf("\n\nThroughput is %d/%d bytes per second\n", CLOCKS_PER_SEC*count*4/(end-start), CLOCKS_PER_SEC*count*4/(end2-start2));

	fclose(out);
//	return 0;
	system("PAUSE");
}    
        