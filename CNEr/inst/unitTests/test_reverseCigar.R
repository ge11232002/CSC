

test_reverseCigar = function(){
  cigar = "16I20M17D"
  checkIdentical(reverseCigar(cigar), "17D20M16I")
}

