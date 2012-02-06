#/bin/bash
patch -p0 < gres_c.patch
patch -p0 < gres_h.patch
patch -p0 < job_scheduler.patch
patch -p0 < makefile.patch
cp -r plugins/sched/lpsched ../src/plugins/sched
cp -r plugins/select/lpconsres ../src/plugins/select
echo "installation completed. run autogen.sh, configure, and then make install!"


