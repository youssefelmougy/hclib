smith_waterman:

```make``` generates three executable files: smith_waterman.baseline.out, smith_waterman.replay.out, smith_waterman.replication.out.

Command-line options:

-tile_width <int>: Width of a single tile (default 32)

-tile_height <int>: Height of a single tile (default 32)

-[no]checksum: Whether or not to use checksums (default checksum)

-inject <filename>: Failure list for failure injection (NOT Supported)

-inject_rate <filename>: Probability of failures

-input1 <filename>: Input file 1

-input2 <filename>: Input file 2

run.sh:

```
./run.sh <WORKLOAD> <APP>
```
