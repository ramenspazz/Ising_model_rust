use std::cmp::Ordering;

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
	match cap.cmp(&0_usize) {
	    Ordering::Equal => {
	        let (tx, rx) = crossbeam_channel::unbounded();
                Self { tx, rx }
	    },
            Ordering::Greater => {
                let (tx, rx) = crossbeam_channel::bounded(cap);
                Self { tx, rx }
            },
            Ordering::Less => panic!("cap must be positive or equal to 0!"),
	}
    }

    pub fn close_channel(self) {
        drop(self.tx);
    }

    pub fn send(&self, msg: T) -> Result<(), crossbeam_channel::SendError<T>> {
        self.tx.send(msg)
    }

    /// Blocks the current thread until a message is received or the channel is empty and disconnected.
    /// If the channel is empty and not disconnected, this call will block until the receive operation can proceed. If the channel is empty and becomes disconnected, this call will wake up and return an error.
    /// If called on a zero-capacity channel, this method will wait for a send operation to appear on the other side of the channel.
    pub fn recv(&self) -> Result<T, crossbeam_channel::RecvError> {
        self.rx.recv()
    }

    ///Attempts to receive a message from the channel without blocking.
    ///This method will either receive a message from the channel immediately or return an error if the channel is empty.
    ///
    ///If called on a zero-capacity channel, this method will receive a message only if there happens to be a send operation on the other side of the channel at the same time.
    pub fn try_recv(&self) -> Result<T, crossbeam_channel::TryRecvError> {
        self.rx.try_recv()
    }

    pub fn recv_timeout(
        &self,
        timeout: std::time::Duration,
    ) -> Result<T, crossbeam_channel::RecvTimeoutError> {
        self.rx.recv_timeout(timeout)
    }

    pub fn is_full(&self) -> bool {
        self.rx.is_full()
    }
}
