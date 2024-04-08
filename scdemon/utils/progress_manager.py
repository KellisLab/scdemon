
class ProgressManager:
    def __init__(self):
        self.pbar = None
    def hook(self, n):
        from tqdm.auto import tqdm
        if self.pbar is None:
            self.pbar = tqdm(total=n)
        elif n < 0:
            self.pbar.close()
            self.pbar = None
        else:
            self.pbar.update(n)

def _interrupt_checker():
    pass
