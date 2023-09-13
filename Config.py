import yaml

class Config(object):
    def __init__(self, path = None):
        self.path = path

        if self.path is not None:
            self.read_file()

    def set_path(self, path):
        self.path = path

    def read_file(self):
        with open(self.path, 'r') as f:
            self.v = yaml.safe_load(f)