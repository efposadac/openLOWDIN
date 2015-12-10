# Color string
def str_green(string):
    return chr(27) + "[1;32m" + string + chr(27) + "[0m"


def str_red(string):
    return chr(27) + "[1;31m" + string + chr(27) + "[0m"
