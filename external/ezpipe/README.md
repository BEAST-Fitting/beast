EZPIPE -- A simple pipeline using coroutines
============================================

When the logic for a function with complex behavior is divided into several
self-contained steps that are themselves functions, these functions are called
*helper functions* or *subroutines*. Subroutines are called by a main function
that is responsible for coordinating the use of several subroutines.

    main-func
        sub-func-1
        sub-func-2
        ...
        sub-funct-n

An approach that is particularly applicable to the task of processing
sequential data, is to use Coroutines. Like a subroutine, a coroutine computes
a single step of a complex computation. However, **when using coroutines**, there
is **no main function** to coordinate results. Instead coroutines themselves link
together to form a **pipeline**.

Coroutines produce products that they send to other coroutines. Processing data
from end to end is a succession of steps or coroutines to finally obtain the
final result.
Coroutines are like a **succesion of asynchrone iterators/generators.**
They cooperate to form a pipeline without any supervising function responsible
for calling them.  This becomes very handy when processing a lot of different
inputs with the same code.
**Coroutines are like workers, which after an initialization, waits for some
job to do** and they can **work in parallel**, while a function works on some
data, the producer can prepare the next iteration.

                     > work1 >
    produce >> work >> work2 >> work >> consume
                     > work3 >


* A Producer is the initial task of a pipeline.  It creates items in a series and uses send(), but not (yield)
* A Filter uses (yield) to consume items and send() to send result to a next step.
* A Consumer uses (yield) to consume items, but does not send.

A deep description is given at `lecture1`_, a simple tutorial at `lecture2`_,
many exmaples are `here`_

.. _lecture1: http://wla.berkeley.edu/~cs61a/fa11/lectures/streams.html
.. _lecture2: http://sdiehl.github.com/coroutine-tutorial/
.. _here: http://www.dabeaz.com/coroutines/

