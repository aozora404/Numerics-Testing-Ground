#include <stdio.h>
#include <stdlib.h>

int main() {

    int x = 0;
    int a = 76;
    int b = 664;
    int c = 296;

    while (a*x % b != c) {
        x ++;
    }
    printf("%d*%d = %d (mod %d)\n", a, x, (a * x) % b, b);

    return 0;
}