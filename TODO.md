* Run rustfmt on all code and make sure code is well formatted. Ignore legacy.
* Write a YAML summary of each prop design, with details on performance RPM e.t.c, as well as section lists. This would default to <propname.yml> in the output directory.
* Make sure polar seeding warm up uses multiple cores if possible.
* The phase after warm up takes quite a while and produces no output. Please use progress bar there.
* composed-camber warm-up should use multiple cores if possible.
* Use a progress bar for the phase after composed-camber-warmup.
