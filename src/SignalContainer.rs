#[derive(Clone, PartialEq)]
pub enum SignalType {
    SigStart,
    SigStop,
}

#[derive(Clone)]
pub struct SignalContainer<T> {
    tx: crossbeam_channel::Sender<T>,
    rx: crossbeam_channel::Receiver<T>,
}

impl<T> SignalContainer<T> {
    /// Creates a channel of bounded capacity.
    ///
    /// This channel has a buffer that can hold at most `cap` messages at a time.
    ///
    /// A special case is zero-capacity channel, which cannot hold any messages. Instead, send and
    /// receive operations must appear at the same time in order to pair up and pass the message over.
    pub fn new(cap: usize) -> Self {
        if cap == 0 {
            let (tx, rx) = crossbeam_channel::unbounded();
            return Self {
                tx,
                rx,
            }
        }
        else if cap > 0 {
            let (tx, rx) = crossbeam_channel::bounded(cap);
            return Self {
                tx,
                rx,
            }
        }
        panic!("Error occured!");
    }
    
    pub fn close_channel(self) { drop(self.tx); }

    pub fn send(&self, msg: T) -> Result<(), crossbeam_channel::SendError<T>> { self.tx.send(msg) }

    pub fn recv(&self) -> Result<T, crossbeam_channel::RecvError> { self.rx.recv() }

    pub fn try_recv(&self) -> Result<T, crossbeam_channel::TryRecvError> { self.rx.try_recv() }

    pub fn is_full(&self) -> bool { self.rx.is_full() }
}