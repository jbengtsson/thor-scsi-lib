'''Observers to be used together with thor_scsi elements

See :class:`Observer`
'''
from .lib import (ObservedState, Observer as _AbstractObserver,
                  ss_vect_tps_to_lists)


class Observer(_AbstractObserver):
    '''Observer example

    Warning:
       Currently it was only tested for ss_vect_tps

    '''
    def __init__(self):
        _AbstractObserver.__init__(self)
        self.name = None
        self.index = None
        self.ps = None

    def __repr__(self):
        cls_name = self.__class__.__name__
        txt = f'{cls_name}(name={self.name}, index={self.index}, ps={self.ps})'
        return txt

    def reset(self):
        self.ps = None

    def view(self, element, ps, observed_state, cnt):
        '''Current view of state at eleement
        Args:
            element:        the element that is observed
            ps:             phase space state
            observed_state: at which state the elemnt is observed
            cnt:            some internal count (e.g. integration step)
                            can be used freely by the element
                            intended to be used with ObservedState.event

        The observed_state and cnt are a bit unusual. These were inspired
        by bluesky's event_document model. Furthermore it can be useful to
        view internal state of an integrator.
        '''
        if observed_state == ObservedState.start:
            self.reset()
            if not self.name:
                self.name = element.name
                self.index = element.index
            return

        elif observed_state == ObservedState.end:
            # Memory management to be reviewed ...
            self.ps = ss_vect_tps_to_lists(ps)

        # Other observed states are not recognised
