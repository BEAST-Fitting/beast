""" Pipeline core
TODO:
* add task group:
    task group is an unmanaged pipeline, or a big task the goal is to allow
    parallel tasks with identical inputs, and improvement to multi-process
    parts this will require a reduce task a the end.

    This must be like a big task Maybe using operator + on Tasks to
    generate a group?
"""

import sys
from .helpers import Task, Logger


class Pipeline(object):
    """
    * a pipeline is a succesion of tasks that are dependent or independent
    * a task is a (group of) function(s) that consume an input and produce an output

    The initial or final tasks are normal tasks, i.e., they both take inputs and could return outputs.
    My advice is to make the first task a checkpoint

    intermediate tasks are like filters: they process data and propagate results.
    """
    def __init__(self, name=None, tasks=(), logger=None, **kwargs):
        """ Constructor
        KEYWORDS
        --------
        name: str
            pipeline name, useful with logs and disk storage
            default is 'pipeline'

        tasks:   seq
            sequence of Task objects

        logger:  buffer/Logger/str
            Logger initialization (see set_logger for help)
        """
        if logger is None:
            logger = sys.stdout
        self.name = name or 'pipeline'
        self.tasks = []
        self.task_names = []
        self._njobs = 0
        for k in tasks:
            self.append(k)
        self.set_logger(logger)

    def get_task(self, task_name):
        """ return pipeline tasks
        INPUTS
        ------
        task_name:   sequence(str/int)
            task identification by name or index
        """
        if hasattr(task_name, '__iter__'):
            return [self.get_task(k) for k in task_name]
        else:
            if type(task_name) == int:
                if task_name > len(self.tasks):
                    raise IndexError('List index {} out of range.'.format(task_name))
                return self.tasks[task_name]
            elif type(task_name) == str:
                if task_name not in self.task_names:
                    raise NameError('Taks {} not found.'.format(task_name))
                for e, k in enumerate(self.task_names):
                    if k == task_name:
                        return self.tasks[e]

    def clear_memoized(self, task_name=None, clearall=False):
        """ Clear memoized values
        KEYWORDS
        --------
        task_name:   sequence(str/int)
            task identification by name or index
            if None, all tasks will receive the update

        clearall:    bool
            if set all tasks will be cleared
        """
        if clearall:
            tlist = self.tasks
        else:
            if task_name is None:
                return
            else:
                tlist = self.get_task(task_name)
                for k in tlist:
                    k.clear_memoized()

    def set_memoized(self, val, task_name=None):
        """ Update task memoization status individually or as a whole

        INPUTS
        ------
        val:  bool/str
            if False, no storage is involved
            if True, storage will be in memory
            if set with a string, PickleShareDB will be invoved (see also Task.memoize)

        KEYWORDS
        --------
        task_name:   sequence(str/int)
            task identification by name or index if None, all tasks will receive the update
        """
        if type(val) == str:
            _val = '{}/{}'.format(self.name, val)
        else:
            _val = val
        if task_name is None:
            tlist = self.tasks
        else:
            tlist = self.get_task(task_name)
        for k in tlist:
            k.memoize(_val)

    def set_logger(self, logger, task_name=None, **kwargs):
        """ Update Logger and propagate the update to the individual tasks
        INPUTS
        ------
        logger:  Logger/str/buffer
            Logger class instanciation parameter

        KEYWORDS
        --------
        task_name:   sequence(str/int)
            task identification by name or index
            if None, all tasks will receive the update

        **kwargs    forwarded to Logger.__init__()
        """
        if task_name is None:
            tlist = self.tasks
        else:
            tlist = self.get_task(task_name)
        for k in tlist:
            k.set_logger(logger, keep_opened=False, **kwargs)

        if issubclass(logger.__class__, Logger):
            self.logger = logger
        else:
            self.logger = Logger(out=logger, **kwargs)

    def __get_next_jobid__(self):
        """ internal
        incrementation of the job counter and name generation for a job
        """
        self._njobs += 1
        return '{}_{}'.format(self.name, self._njobs)

    def append(self, t):
        """ Append a task to the pipeline
        INPUTS
        ------
        t:   Task
            task to add
        """
        if not issubclass(t.__class__, Task):
            raise ValueError('Expecting a Task instance, got {}'.format(type(t)))
        self.tasks.append(t)
        self.task_names.append(t.name)

    def insert(self, idx, t):
        """ Insert a task to the pipeline at a given place
        INPUTS
        ------
        idx: int
            place to insert the task in the line

        t:  Task
            task to add
        """
        if not issubclass(t.__class__, Task):
            raise ValueError('Expecting a Task instance, got {}'.format(type(t)))
        self.tasks.insert(idx, t)
        self.task_names.insert(idx, t.name)

    def __call__(self, *args, **kwargs):
        """ Starts the pipeline from the first task with the 'arg' arguments"""
        job_id = self.__get_next_jobid__()
        return self.start_job_from(0, job_id, val=args)

    def __repr__(self):
        """ Str representation of the pipeline """
        txt = 'Pipeline: {}, {}\n'.format(self.name, object.__repr__(self))
        for e, k in enumerate(self.tasks):
            memo = 'M' if k.memoized is not False else '-'
            memo += 'L' if k.logger.out is not None else '-'
            txt += '   {} [{}] {}\n'.format(e, memo, k.name)
        return txt

    def get_job_status(self, job_id):
        """ Return the status of a given job through each task
        It will return True each time a job went through a task and left a memoized value
        This is mainly used to resume the jobs or check if the job went through.

        INPUTS
        ------
        job_id:  str
            job identifier
        """
        _job_id = '{}_{}'.format(self.name, job_id)
        return [ k.get_job_status(_job_id, strict=False) for k in self.tasks ]

    def resume_job(self, job_id, verbose=True):
        """ resume a job if possible. It will check a job status through the
        line and figure where the last checkpoint was to restart from there

        INPUTS
        ------
        job_id:  str
            job identifier

        KEYWORDS
        --------
        verbose: bool
            if set, prints out the status
        """

        stat = self.get_job_status(job_id)
        if sum(stat) == 0:
            raise LookupError('No available checkpoint to resume Job "{}"'.format(job_id))
        else:
            # find the last available restart
            for k in range(len(self.tasks))[::-1]:
                if stat[k]:
                    break
            if k == len(self.tasks) - 1:
                if verbose:
                    print('Job {} already completed'.format(job_id))
                return
            if verbose:
                print('Restarting Job "{}" from Task "{}"'.format(job_id, self.task_names[k]))
                # call the last memoized task and restart from there
                self.start_job_from(k, job_id, val=None)

    def start_job_from(self, task_name, job_id, val=None):
        """ Start the pipeline from any Task
        INPUTS
        ------
        task_name:   int/str
            task indice or name

        job_id:      str
            job identification name

        KEYWORDS
        --------
        val:  object
            value to pipe to the first task of the sequence
        """
        if type(task_name) == int:
            if task_name > len(self.tasks):
                raise IndexError('List index {} out of range.'.format(task_name))
            tlist = self.tasks[task_name:]
        elif type(task_name) == str:
            if task_name not in self.task_names:
                raise NameError('Taks {} not found.'.format(task_name))
            for e, k in enumerate(self.task_names):
                if k == task_name:
                    break
            tlist = self.tasks[e:]

        _val = (job_id, val)
        for tk in tlist:
            _val = tk.__ror__(_val)
        return _val

    def __ror__(self, other):
        """ Make a pipeline pipeable, a pipeline could be considered a one task
        INPUTS
        ------
        other:   object
            value to send to the first task

        OUTPUTS
        -------
        out:  object
            output from the last task of the pipeline
        """
        return self(other)
