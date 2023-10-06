# openMP-labwork
## Notes:
- include omp.h
- multi-threading, shared address model - threads communicate by sharing variables
- unintended sharing of data causes *race conditions*(program's outcome changes as the threads are scheduled differently)
- use synchronization to protect data conflicts
- synchronization - expensive - change how data is accesses to minimize the need for it

### Fork-join parallelism
- master thread spawns a team of threads as needed
- nested parallel region joined by sequential parts via master thread
- imagine a series of parallel resistor networks

### Thread creation
```
double A[100];
omp_set_num_threads(4); //runtime fxn to request #N threads
#pragma omp parallel
{
int ID = omp_get_thread_num(); //returns thread ID
foobar(ID,A); //each thread call this function for ID=0,1,2,3
}
std::cout << "All done";
```
### Synchronisation
`#pragma omp sync_type`

*High Level synchronisation*
- critical: mutual exclusion - only 1 thread at a time can enter a critical region
- atomic
- barrier: each thread waits until all threads arrive
- ordered

*Low level sync*
- flush
- locks (both simple and nested)

TO-DO: update the doc
