using CsvHelper;
using Mondrian.Metadata;
using Mondrian.Models.Input;
using Mondrian.Models.Source;
using System.Globalization;
using System.Text.RegularExpressions;

Console.WriteLine("Opening Source Csv..");
using (var reader = new StreamReader("metadata_scy-263.csv"))
using (var csv = new CsvReader(reader, CultureInfo.InvariantCulture))
{
  Console.WriteLine("Parsing Csv..");
  var prefix = "/rrdevcfa2047244445/inputs/scy-263/dataset/SCY-263/";
  var records = csv.GetRecords<LeapCsv>();
  var cellRecords = new Dictionary<CellRecord, InputCell>();
  foreach (var record in records)
  {
    var pathParts = record.file_paths.Split('/');
    var pathTrunc = new ArraySegment<string>(pathParts, 6, pathParts.Length - 6);
    var azurePath = prefix + string.Join('/', pathTrunc);
    string filename = Path.GetFileName(record.file_paths);
    string pattern = @"(.*C\d{1,2}_)(S_\d).(R{1,2}_\d).fq.gz";
    Regex regex = new(pattern);
    Match match = regex.Match(filename);
    CellRecord cellRecord = new(record.cell, record.flow_cells, match.Groups[2].Value);
    if (cellRecords.ContainsKey(cellRecord))
    {
      switch (match.Groups[3].Value)
      {
        case "R_1":
          cellRecords[cellRecord].Lane.fastq1 = azurePath;
          break;
        case "R_2":
          cellRecords[cellRecord].Lane.fastq2 = azurePath;
          break;
        default:
          break;
      }
    }
    else
    {
      if (match.Success)
      {
        var Cell = new Mondrian.Models.Aggregate.Cell()
        {
          cell_id = record.cell,
          column = record.column,
          condition = "A",
          is_control = record.is_control == "TRUE",
          library_id = record.library_id,
          primer_i5 = record.index_i5_list,
          primer_i7 = record.index_i7_list,
          row = record.row,
          sample_id = record.sample_id,
          sample_type = "A",
          lanes = new List<Mondrian.Models.Input.Lane>()
        };
        switch (match.Groups[3].Value)
        {
          case "R_1":
            cellRecords.Add(cellRecord,
            new InputCell()
            {
              CellId = record.cell,
              Cell = (Mondrian.Models.Aggregate.Cell)Cell,
              Lane = new Mondrian.Models.Aggregate.Lane()
              {
                fastq1 = azurePath,
                fastq2 = "",
                flowcell_id = record.flow_cells,
                lane_id = cellRecord.Lane
              }
            }
            );
            break;
          case "R_2":
            cellRecords.Add(cellRecord,
            new InputCell()
            {
              CellId = record.cell,
              Cell = (Mondrian.Models.Aggregate.Cell)Cell,
              Lane = new Mondrian.Models.Aggregate.Lane()
              {
                fastq1 = azurePath,
                fastq2 = "",
                flowcell_id = record.flow_cells,
                lane_id = cellRecord.Lane
              }
            });
            break;
          default:
            break;
        }
      }
    }
  }
  
  List<InputCell> InputCells = cellRecords.Values.ToList();
  // Generate metadata.yaml and write to file
  Console.WriteLine("Writing metadata.yaml");
  Metadata metadata = new(InputCells);
  using (StreamWriter outputFile = new("scy263-metadata.yaml"))
  {
    await outputFile.WriteAsync(metadata.Yaml);
  }
  // Generate inputs.json and write to file
  AlignmentWorkflow alignmentWorkflow = new()
  {
    docker_image = "quay.io/mondrianscwgs/alignment:v0.0.84",
    metadata_yaml = "/rrdevcfa2047244445/inputs/scy-263/scy263-metadata.yaml",
    reference = new ReferenceGenome()
    {
      genome_name = "human",
      reference = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/human/GRCh37-lite.fa",
      reference_fa_fai = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/human/GRCh37-lite.fa.fai",
      reference_fa_amb = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/human/GRCh37-lite.fa.amb",
      reference_fa_ann = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/human/GRCh37-lite.fa.ann",
      reference_fa_bwt = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/human/GRCh37-lite.fa.bwt",
      reference_fa_pac = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/human/GRCh37-lite.fa.pac",
      reference_fa_sa = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/human/GRCh37-lite.fa.sa",
    },
    supplimentary_references = new List<ReferenceGenome>
  {
    new ReferenceGenome()
    {
      genome_name = "mouse",
      reference = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/mouse/mm10_build38_mouse.fasta",
      reference_fa_fai = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/mouse/mm10_build38_mouse.fasta.fai",
      reference_fa_amb = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/mouse/mm10_build38_mouse.fasta.amb",
      reference_fa_ann = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/mouse/mm10_build38_mouse.fasta.ann",
      reference_fa_bwt = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/mouse/mm10_build38_mouse.fasta.bwt",
      reference_fa_pac = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/mouse/mm10_build38_mouse.fasta.pac",
      reference_fa_sa = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/mouse/mm10_build38_mouse.fasta.sa"
    },
    new ReferenceGenome()
    {
      genome_name = "salmon",
      reference = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/salmon/GCF_002021735.1_Okis_V1_genomic.fna",
      reference_fa_fai = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/salmon/GCF_002021735.1_Okis_V1_genomic.fna.fai",
      reference_fa_amb = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/salmon/GCF_002021735.1_Okis_V1_genomic.fna.amb",
      reference_fa_ann = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/salmon/GCF_002021735.1_Okis_V1_genomic.fna.ann",
      reference_fa_bwt = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/salmon/GCF_002021735.1_Okis_V1_genomic.fna.bwt",
      reference_fa_pac = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/salmon/GCF_002021735.1_Okis_V1_genomic.fna.pac",
      reference_fa_sa = "/rrdevcfa2047244445/datasets/reference/mondrian-ref-GRCh37/salmon/GCF_002021735.1_Okis_V1_genomic.fna.sa"
    }
  },
    fastq_files = new List<Mondrian.Models.Input.Cell>()
  };
  Console.WriteLine("Writing inputs.json");
  Inputs inputs = new(InputCells, alignmentWorkflow);
  using (StreamWriter outputFile = new("scy263-inputs.json"))
  {
    await outputFile.WriteAsync(inputs.Json);
  }
}
