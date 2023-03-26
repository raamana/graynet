
from sys import version_info

if version_info.major > 2:
    from graynet import run_workflow
else:
    raise NotImplementedError('Python 3 or higher is required to run graynet.'
                              'Please upgrade.')

def main():
    "Entry point."

    run_workflow.cli_run()

if __name__ == '__main__':
    main()
