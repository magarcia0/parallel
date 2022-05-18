#include<stdio.h>
#include<math.h>



int main() {
  int soln =0;
  for(int i = 0; i <25; i++){
    soln = i % 5;
    printf("5 mod %d = %d\n", i, soln);
  }

  return 0;
}
