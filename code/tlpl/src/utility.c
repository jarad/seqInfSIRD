
/* Utility function to determine if any element in a vector is negative */
int anyNegative(int n, int *v)                     // both should be const
{
    int i;
    for (i=0; i<n; i++) 
    {
        if (v[i]<0) return 1;
    }
    return 0;
}


/* Integer vector copy */
void copy_int(int n, int *from, int *to) 
{
    int i;
    for (i=0; i<n; i++) to[i]=from[i];
}
