"""Wait at most 30 minutes for our older queued file; retry the final upload once."""
import ftplib
import pathlib
import subprocess
import time

deadline = time.monotonic() + 1800
target = '/R-devel/mixeff_0.2.0.tar.gz'
original = '213 20261026194051'
print('Waiting for the earlier win-builder R-devel upload to leave the queue.', flush=True)
while time.monotonic() < deadline:
    try:
        with ftplib.FTP('win-builder.r-project.org', timeout=20) as ftp:
            ftp.login()
            try:
                stamp = ftp.sendcmd('MDTM ' + target)
            except ftplib.error_perm as error:
                if not str(error).startswith('550'):
                    raise
                print('Earlier queue entry no longer exposes an MDTM; trying one upload.', flush=True)
                with open('win-R-devel-retry.log', 'w') as log:
                    result = subprocess.run(
                        ['Rscript', 'submit-builder.R', 'win-R-devel'],
                        stdout=log, stderr=subprocess.STDOUT, timeout=180,
                    )
                print('Final upload exit status:', result.returncode, flush=True)
                raise SystemExit(result.returncode)
            if stamp != original:
                print('Queue entry changed unexpectedly; no upload attempted:', stamp, flush=True)
                raise SystemExit(2)
    except (OSError, ftplib.error_temp) as error:
        print('Temporary metadata check failure:', str(error), flush=True)
    time.sleep(min(60, max(0, deadline - time.monotonic())))
print('Wait budget exhausted; no retry attempted.', flush=True)
raise SystemExit(3)
