#!/usr/bin/perl

my $file1 = $ARGV[0];
my $col1 = $ARGV[1];
my $file2 = $ARGV[2];
my $col2 = $ARGV[3];
my $neg = $ARGV[4];

my %present = ();

open (F1, $file1);
while(<F1>){
    chomp;
    @fields = split;
    $present{ $fields[$col1-1] } = 1;
}
close F1;

open (F2, $file2);
while (<F2>){
    chomp;
    @fields = split;
    if ( ( $neg eq "-v" && ! exists $present{$fields[$col2-1]} ) || ( $neg eq "" && exists $present{$fields[$col2-1]} ) ) { print "$_\n"; }
}
close F2;


