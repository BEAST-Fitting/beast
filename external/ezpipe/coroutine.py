"""

* a pipeline is a succesion of tasks that are dependent or independent
* a task is a (group of) function(s) that consume an input and produce an output

Based on Gist: https://gist.github.com/sampsyo/502867

Simple but robust implementation of generator/coroutine-based pipelines in
Python. The pipelines may be run either sequentially (single-threaded) or in
parallel (one thread per pipeline stage).

This implementation supports pipeline bubbles (indications that the
processing for a certain item should abort). To use them, yield the
BUBBLE constant from any stage coroutine except the last.

In the parallel case, the implementation transparently handles thread
shutdown when the processing is complete and when a stage raises an
exception.
"""
from queue import Queue
from threading import Thread, Lock

BUBBLE = '__PIPELINE_BUBBLE__'
POISON = '__PIPELINE_POISON__'

DEFAULT_QUEUE_SIZE = 100

__all__ = ['Pipeline']


class PipelineError(object):
    """An indication that an exception occurred in the pipeline. The
    object is passed through the pipeline to shut down all threads
    before it is raised again in the main thread.
    """
    def __init__(self, exc):
        self.exc = exc


class FirstPipelineThread(Thread):
    """The thread running the first stage in a parallel pipeline setup.
    The coroutine should just be a generator.
    """
    def __init__(self, coro, out_queue):
        super(FirstPipelineThread, self).__init__()
        self.coro = coro
        self.out_queue = out_queue

        self.abort_lock = Lock()
        self.abort_flag = False

    def run(self):
        while True:
            # Time to abort?
            with self.abort_lock:
                if self.abort_flag:
                    break

            # Get the value from the generator.
            try:
                msg = next(self.coro)
            except StopIteration:
                break
            except Exception as exc:
                self.out_queue.put(PipelineError(exc))
                return

            # Send it to the next stage.
            self.out_queue.put(msg)
            if msg is BUBBLE:
                continue

        # Generator finished; shut down the pipeline.
        self.out_queue.put(POISON)

    def abort(self):
        """Shut down the pipeline by canceling this thread and
        poisoning out_channel.
        """
        with self.abort_lock:
            self.abort_flag = True


class MiddlePipelineThread(Thread):
    """A thread running any stage in the pipeline except the first or
    last.
    """
    def __init__(self, coro, in_queue, out_queue):
        super(MiddlePipelineThread, self).__init__()
        self.coro = coro
        self.in_queue = in_queue
        self.out_queue = out_queue

    def run(self):
        # Prime the coroutine.
        next(self.coro)

        while True:
            # Get the message from the previous stage.
            msg = self.in_queue.get()
            if msg is POISON:
                break
            elif isinstance(msg, PipelineError):
                self.out_queue.put(msg)
                return

            # Invoke the current stage.
            try:
                out = self.coro.send(msg)
            except Exception as exc:
                self.out_queue.put(PipelineError(exc))
                return

            # Send message to next stage.
            if out is BUBBLE:
                continue
            self.out_queue.put(out)

        # Pipeline is shutting down normally.
        self.out_queue.put(POISON)


class LastPipelineThread(Thread):
    """A thread running the last stage in a pipeline. The coroutine
    should yield nothing.
    """
    def __init__(self, coro, in_queue):
        super(LastPipelineThread, self).__init__()
        self.coro = coro
        self.in_queue = in_queue

    def run(self):
        # Prime the coroutine.
        next(self.coro)

        while True:
            # Get the message from the previous stage.
            msg = self.in_queue.get()
            if msg is POISON:
                break
            elif isinstance(msg, PipelineError):
                self.exc = msg.exc
                return

            # Send to consumer.
            try:
                self.coro.send(msg)
            except Exception as exc:
                self.exc = exc
                return

        # No exception raised in pipeline.
        self.exc = None


class Pipeline(object):
    """Represents a staged pattern of work. Each stage in the pipeline
    is a coroutine that receives messages from the previous stage and
    yields messages to be sent to the next stage.
    """
    def __init__(self, stages):
        """Makes a new pipeline from a list of coroutines. There must
        be at least two stages.
        """
        if len(stages) < 2:
            raise ValueError('pipeline must have at least two stages')
        self.stages = stages

    def run_sequential(self):
        """Run the pipeline sequentially in the current thread. The
        stages are run one after the other.
        """
        # "Prime" the coroutines.
        for coro in self.stages[1:]:
            next(coro)

        # Begin the pipeline.
        for msg in self.stages[0]:
            for stage in self.stages[1:]:
                msg = stage.send(msg)
                if msg is BUBBLE:
                    # Don't continue to the next stage.
                    break

    def run_parallel(self, queue_size=DEFAULT_QUEUE_SIZE):
        """Run the pipeline in parallel using one thread per stage. The
        messages between the stages are stored in queues of the given
        size.
        """
        queues = [Queue(queue_size) for i in range(len(self.stages) - 1)]
        threads = [FirstPipelineThread(self.stages[0], queues[0])]
        for i in range(1, len(self.stages) - 1):
            threads.append(MiddlePipelineThread(
                self.stages[i], queues[i - 1], queues[i]
            ))
        threads.append(LastPipelineThread(self.stages[-1], queues[-1]))

        # Start threads.
        for thread in threads:
            thread.start()

        # Wait for termination.
        try:
            for thread in threads:
                thread.join()
        except:
            # Shut down the pipeline by telling the first thread to
            # poison its channel.
            threads[0].abort()
            raise

        # Was there an exception?
        exc = threads[-1].exc
        if exc:
            raise exc


def xxtest_case():
    """Smoke test."""
    import time

    def produce(pid=0):
        for i in range(5):
            print(pid, 'generating', i)
            time.sleep(0.1)
            yield i

    def work(pid):
        num = yield
        while True:
            print(pid, 'waiting...')
            time.sleep(2)
            num = yield num * 2
            print(pid, 'sent', num)

    def consume(pid=-1):
        while True:
            print(pid, 'waiting...')
            num = yield
            time.sleep(1)
            print(pid, 'received', num)

    ts_start = time.time()
    Pipeline([produce(0), work(1), work(2), consume(3)]).run_sequential()
    ts_middle = time.time()
    Pipeline([produce(0), work(1), work(2), consume(3)]).run_parallel()
    ts_end = time.time()

    print('Sequential time:', ts_middle - ts_start)
    print('Parallel time:', ts_end - ts_middle)


if __name__ == '__main__':
    test_case()
