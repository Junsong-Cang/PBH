#include <stdio.h>
#include<math.h>

int main()
{
    int i = 0;

    while (1)
    {
        printf("---%d\n", i);

        if (i > 100)
        {
            break; // Break out of the while loop when i is equal to 5
        }
        printf("%d\n", i);
        i++;
    }

    printf("Loop ended\n");

    double xx = fmin(1.2, 3.1);
    printf("%f\n", xx);

    return 0;
}
