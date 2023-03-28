from sys import version_info

if version_info.major > 2:
    from graynet import cli_run
else:
    raise NotImplementedError('Python 3 or higher is required to run graynet.'
                              'Please upgrade.')
del version_info


def main():
    "Entry point."

    cli_run()


if __name__ == '__main__':
    main()
