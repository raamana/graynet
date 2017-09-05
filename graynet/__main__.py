
from sys import version_info

if version_info.major==2 and version_info.minor==7 and version_info.micro==13:
    import graynet
elif version_info.major > 2:
    from graynet import graynet
else:
    raise NotImplementedError('hiwenet supports only 2.7.13 or 3+. Upgrade to Python 3+ is recommended.')

def main():
    "Entry point."

    graynet.cli_run()

if __name__ == '__main__':
    main()
