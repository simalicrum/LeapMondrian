using CsvHelper;
using Microsoft.Extensions.Configuration;
using Mondrian.Metadata;
using Mondrian.Models.Input;
using Mondrian.Models.Source;
using System.Globalization;
using System.Text.RegularExpressions;

var configuration = new ConfigurationBuilder()
    .SetBasePath(Directory.GetCurrentDirectory())
    .AddJsonFile("appsettings.json")
    .Build();


var inputFilePath = configuration["InputFilePath"];
var outputMetadataYaml = configuration["OutputMetadataYaml"];
var outputInputsJson = configuration["OutputInputsJson"];
var dataStoragePrefix = configuration["DataStoragePrefix"];
var dockerImage = configuration["DockerImage"];
var metadataYamlPath = configuration["MetadataYamlPath"];
var referenceGenomes = configuration["ReferenceGenomes"];

Console.WriteLine("Opening Source Csv..");

using (var reader = new StreamReader(inputFilePath))
using (var csv = new CsvReader(reader, CultureInfo.InvariantCulture))
{
  Console.WriteLine("Parsing Csv..");
  var prefix = dataStoragePrefix;
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
                fastq1 = "",
                fastq2 = azurePath,
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
  using (StreamWriter outputFile = new(outputMetadataYaml))
  {
    await outputFile.WriteAsync(metadata.Yaml);
  }
  // Generate inputs.json and write to file
  AlignmentWorkflow alignmentWorkflow = new()
  {
    docker_image = dockerImage,
    metadata_yaml = metadataYamlPath,
    reference = new ReferenceGenome()
    {
      genome_name = "human",
      reference = $"{referenceGenomes}/human/GRCh37-lite.fa",
      reference_fa_fai = $"{referenceGenomes}/human/GRCh37-lite.fa.fai",
      reference_fa_amb = $"{referenceGenomes}/human/GRCh37-lite.fa.amb",
      reference_fa_ann = $"{referenceGenomes}/human/GRCh37-lite.fa.ann",
      reference_fa_bwt = $"{referenceGenomes}/human/GRCh37-lite.fa.bwt",
      reference_fa_pac = $"{referenceGenomes}/human/GRCh37-lite.fa.pac",
      reference_fa_sa = $"{referenceGenomes}/human/GRCh37-lite.fa.sa",
    },
    supplimentary_references = new List<ReferenceGenome>
  {
    new ReferenceGenome()
    {
      genome_name = "mouse",
      reference = $"{referenceGenomes}/mouse/mm10_build38_mouse.fasta",
      reference_fa_fai = $"{referenceGenomes}/mouse/mm10_build38_mouse.fasta.fai",
      reference_fa_amb = $"{referenceGenomes}/mouse/mm10_build38_mouse.fasta.amb",
      reference_fa_ann = $"{referenceGenomes}/mouse/mm10_build38_mouse.fasta.ann",
      reference_fa_bwt = $"{referenceGenomes}/mouse/mm10_build38_mouse.fasta.bwt",
      reference_fa_pac = $"{referenceGenomes}/mouse/mm10_build38_mouse.fasta.pac",
      reference_fa_sa = $"{referenceGenomes}/mouse/mm10_build38_mouse.fasta.sa"
    },
    new ReferenceGenome()
    {
      genome_name = "salmon",
      reference = $"{referenceGenomes}/salmon/GCF_002021735.1_Okis_V1_genomic.fna",
      reference_fa_fai = $"{referenceGenomes}/salmon/GCF_002021735.1_Okis_V1_genomic.fna.fai",
      reference_fa_amb = $"{referenceGenomes}/salmon/GCF_002021735.1_Okis_V1_genomic.fna.amb",
      reference_fa_ann = $"{referenceGenomes}/salmon/GCF_002021735.1_Okis_V1_genomic.fna.ann",
      reference_fa_bwt = $"{referenceGenomes}/salmon/GCF_002021735.1_Okis_V1_genomic.fna.bwt",
      reference_fa_pac = $"{referenceGenomes}/salmon/GCF_002021735.1_Okis_V1_genomic.fna.pac",
      reference_fa_sa = $"{referenceGenomes}/salmon/GCF_002021735.1_Okis_V1_genomic.fna.sa"
    }
  },
    fastq_files = new List<Mondrian.Models.Input.Cell>()
  };
  Console.WriteLine("Writing inputs.json");
  Inputs inputs = new(InputCells, alignmentWorkflow);
  using (StreamWriter outputFile = new(outputInputsJson))
  {
    await outputFile.WriteAsync(inputs.Json);
  }
}
