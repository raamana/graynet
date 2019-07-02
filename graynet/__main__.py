
from sys import version_info

if version_info.major==2 and version_info.minor==7:
    import run_workflow
elif version_info.major > 2:
    from graynet import run_workflow
else:
    raise NotImplementedError('graynet supports only 2.7.13 or 3+. '
                              'Upgrade to Python 3+ is recommended.')

def main():
    "Entry point."

    run_workflow.cli_run()

if __name__ == '__main__':
    main()
