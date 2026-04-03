function C = elasticityMatrix(beam)

C = diag([beam.GA1, beam.GA2, beam.EA, beam.EI1, beam.EI2, beam.GIt]);