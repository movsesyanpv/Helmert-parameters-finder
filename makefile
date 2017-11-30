run: comp
		./mod
comp: mod.c
		gcc mod.c -lm -llapacke -lproj -o mod
