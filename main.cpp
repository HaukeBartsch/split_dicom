
/*
 ./split_dicom -d data -o /tmp/
*/


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <fstream>
#include <filesystem>

// Some gdcm libraries
#include "gdcmAttribute.h"
#include "gdcmDefs.h"
#include "gdcmGlobal.h"
#include "gdcmImage.h"
#include "gdcmImageChangeTransferSyntax.h"
#include "gdcmImageHelper.h"
#include "gdcmImageReader.h"
#include "gdcmImageWriter.h"
#include "gdcmPhotometricInterpretation.h"
#include "gdcmSystem.h"
#include <gdcmUIDGenerator.h>

// boost libraries
#include "boost/date_time.hpp"
#include "boost/filesystem.hpp"
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "optionparser.h"

struct Arg : public option::Arg {
  static option::ArgStatus Required(const option::Option &option, bool) {
    return option.arg == 0 ? option::ARG_ILLEGAL : option::ARG_OK;
  }
  static option::ArgStatus Empty(const option::Option &option, bool) {
    return (option.arg == 0 || option.arg[0] == 0) ? option::ARG_OK : option::ARG_IGNORE;
  }
};

enum optionIndex { UNKNOWN, INPUT, OUTPUT, HELP };
const option::Descriptor usage[] = {
    {UNKNOWN, 0, "", "", option::Arg::None,
     "Split all input volumes.\n"
     "USAGE: split_volume [options]\n\n"
     "Options:"},
    {HELP, 0, "", "help", Arg::None, "  --help  \tPrint this help message."},
    {INPUT, 0, "i", "input", Arg::Required, "  --input, -i  \tDirectory with DICOM files."},
    {OUTPUT, 0, "o", "output", Arg::Required, "  --output, -o  \tOutput directory."},
    {UNKNOWN, 0, "", "", Arg::None,
     "\nExamples:\n"
     "  ./split_volume -i data/LIDC-IDRI-0009 -o /tmp/bla\n"
     "  ./split_volume --help\n"},
    {0, 0, 0, 0, 0, 0}};

// get a list of all the DICOM files we should use
std::vector<std::string> listFilesSTD(const std::string &path) {
  std::vector<std::string> files;
  std::string extension;

  for (boost::filesystem::recursive_directory_iterator end, dir(path); dir != end; ++dir) {
    // std::cout << *dir << "\n";  // full path
    if (is_regular_file(dir->path())) {
      // reading zip and tar files might take a lot of time.. filter out here
      extension = boost::filesystem::extension(dir->path());
      if (extension == ".tar" || extension == ".gz" || extension == ".zip" || extension == ".tgz" ||
          extension == ".bz2")
        continue;
      files.push_back(dir->path().c_str());
      if ((files.size() % 200) == 0) {
        fprintf(stdout, "[reading files %05lu]\r", files.size());
      }
    }
  }
  fprintf(stdout, "[reading files %'5lu done]\r", files.size());

  return files;
}

int main(int argc, char **argv) {
  setlocale(LC_NUMERIC, "en_US.utf-8");
  // std::srand(std::time(0));
  // std::srand((unsigned)time(NULL) * getpid());

  argc -= (argc > 0);
  argv += (argc > 0); // skip program name argv[0] if present

  option::Stats stats(usage, argc, argv);
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer(stats.buffer_max);
  option::Parser parse(usage, argc, argv, &options[0], &buffer[0]);

  if (parse.error())
    return 1;

  if (options[HELP] || argc == 0) {
    option::printUsage(std::cout, usage);
    return 0;
  }

  std::string input = "";
  std::string output = "";    // directory path
  unsigned int random_seed = (unsigned)time(NULL) * getpid();

  for (option::Option *opt = options[UNKNOWN]; opt; opt = opt->next())
    std::cout << "Unknown option: " << std::string(opt->name, opt->namelen) << "\n";

  for (int i = 0; i < parse.optionsCount(); ++i) {
    option::Option &opt = buffer[i];
    // fprintf(stdout, "Argument #%d is ", i);
    switch (opt.index()) {
    case HELP:
      // not possible, because handled further above and exits the program
    case INPUT:
      if (opt.arg) {
        // fprintf(stdout, "--input '%s'\n", opt.arg);
        input = opt.arg;
      } else {
        fprintf(stdout, "--input needs a directory specified\n");
        exit(-1);
      }
      break;
    case OUTPUT:
      if (opt.arg) {
        // fprintf(stdout, "--output '%s'\n", opt.arg);
        output = opt.arg;
      } else {
        fprintf(stdout, "--output needs a directory specified\n");
        exit(-1);
      }
      break;
    case UNKNOWN:
      // not possible because Arg::Unknown returns ARG_ILLEGAL
      // which aborts the parse with an error
      break;
    }
  }
  std::srand(random_seed);


    //
    // lets read in a random DICOM image (check that it is a DICOM image...)
    //
    bool foundDICOM = false;
    int pickImageIdx;
    int maxIter = 10;
    while (maxIter > 0) {
      pickImageIdx = std::rand() % files.size();
      // to check for DICOM lets use the normal gdcm::Read and CanRead
      gdcm::Reader r;
      r.SetFileName(files[pickImageIdx].c_str());
      if (r.CanRead() != true) {
        maxIter--;
        fprintf(stdout, "read failed [%d] for \"%s\"... try next\n", maxIter,
                files[pickImageIdx].c_str());
        continue; // try again
      }

      reader.SetFileName(files[pickImageIdx].c_str());
      fprintf(stdout, "Try to read \"%s\"\n", files[pickImageIdx].c_str());
      if (!reader.Read()) {
        maxIter--;
        fprintf(stdout, "read failed... try next\n");
        continue; // try again
      }
      foundDICOM = true;
      break;
    }
    if (!foundDICOM) {
      fprintf(stderr, "Error: I tried to find a DICOM image I could read. I could not... ");
      exit(-1);
    } else {
      fprintf(stdout, "  DICOM file: \"%s\"\n", files[pickImageIdx].c_str());
      // bbox.insert(std::pair<std::string, std::string>("filename", files[pickImageIdx].c_str()));
    }
    const gdcm::Image &iimage = reader.GetImage();
    // we should probably change the transfer syntax here, uncompress JPEG2000 so that
    // we can write out the image pair uncompressed...
    gdcm::ImageChangeTransferSyntax change;
    // change.SetTransferSyntax(gdcm::TransferSyntax::JPEG2000Lossless);
    // change.SetTransferSyntax(gdcm::TransferSyntax::JPEGLosslessProcess14_1);
    change.SetTransferSyntax(gdcm::TransferSyntax::ImplicitVRLittleEndian);
    change.SetInput(iimage);
    bool b = change.Change();
    if (!b) {
      std::cerr << "Could not change the Transfer Syntax" << std::endl;
      return 1;
    }
    // unsigned int n_dim = image.GetNumberOfDimensions();
    // const unsigned int *dims = image.GetDimensions();
    // Origin
    // const double *origin = image.GetOrigin();

    gdcm::File &file = reader.GetFile();
    gdcm::DataSet &ds = file.GetDataSet();
    const gdcm::Image &change_image = change.GetOutput();
    const gdcm::PhotometricInterpretation &pi = change_image.GetPhotometricInterpretation();
    // We need to check what the phometric interpretation is. We cannot handle all of them.
    if (pi != gdcm::PhotometricInterpretation::MONOCHROME2) {
      fprintf(stderr, "Warning: We only understand MONOCHROME2, skip this image.\n");
      // we could change the interpreation as well
      /*        gdcm::ImageChangePhotometricInterpretation icpi;
         icpi.SetInput(image);
         icpi.SetPhotometricInterpretation(gdcm::PhotometricInterpretation::MONOCHROME2);
         if (!icpi.Change())
         {
           itkExceptionMacro(<< "Failed to change to Photometric Interpretation");
         }
         itkWarningMacro(<< "Converting from MONOCHROME1 to MONOCHROME2 may impact the meaning of
         DICOM attributes related " "to pixel values."); image = icpi.GetOutput(); */

      continue;
    }
    unsigned short xmax;
    unsigned short ymax;

    // image dimensions
    std::vector<unsigned int> extent = gdcm::ImageHelper::GetDimensionsValue(reader.GetFile());

    xmax = (unsigned short)extent[0];
    ymax = (unsigned short)extent[1];

    gdcm::PixelFormat pf = change_image.GetPixelFormat();
    // pf could be gdcm::PixelFormat::INT8, or INT16, etc...
    // We need to change our cast further down on the pixel data.

    unsigned short pixelsize = pf.GetPixelSize();
    fprintf(stdout, "  pixelsize of input is: %d\n", pixelsize);
    unsigned long len = change_image.GetBufferLength();
    // fprintf(stdout, "Found buffer of length: %ld\n", len);
    char *buffer = new char[len];

    change_image.GetBuffer(buffer);

    // what is max and min here? (don't change by the overlay)
    float current_image_min_value, current_image_max_value;
    {
      signed short *bvals = (signed short *)buffer;
      current_image_min_value = (float)bvals[0];
      current_image_max_value = (float)bvals[0]; 
      for (int p = 1; p < xmax*ymax; p++) {
        if (current_image_min_value > bvals[p]) current_image_min_value = bvals[p];
        if (current_image_max_value < bvals[p]) current_image_max_value = bvals[p];  
      }
    }

    //
    // write the image pair
    //
    dn = output;
    // struct stat buf;
    if (!(stat(dn.c_str(), &buf) == 0)) {
      mkdir(dn.c_str(), 0777);
    }

    dn = output + "/with_text/";
    // struct stat buf;
    if (!(stat(dn.c_str(), &buf) == 0)) {
      mkdir(dn.c_str(), 0777);
    }
    dn = output + "/without_text/";
    if (!(stat(dn.c_str(), &buf) == 0)) {
      mkdir(dn.c_str(), 0777);
    }

    char outputfilename[1024];

      gdcm::ImageWriter writer2;
      writer2.SetImage(change.GetOutput());
      writer2.SetFile(reader.GetFile());
      snprintf(outputfilename, 1024 - 1, "%s/without_text/%08d.dcm", output.c_str(), i);
      writer2.SetFileName(outputfilename);
      if (!writer2.Write()) {
        return 1;
      }


    // now we can add the bitmap to the original data and write again
    // change_image.SetBuffer(buffer);
    // gdcm::DataElement pixeldata = change_image.GetDataElement();
    gdcm::DataElement pixeldata(gdcm::Tag(0x7fe0, 0x0010));
    pixeldata.SetByteValue(buffer, len);

    gdcm::SmartPointer<gdcm::Image> im = new gdcm::Image;
    im->SetNumberOfDimensions(2);
    im->SetDimension(0, xmax);
    im->SetDimension(1, ymax);
    im->SetPhotometricInterpretation(change_image.GetPhotometricInterpretation());
    im->GetPixelFormat().SetSamplesPerPixel(1);
    im->GetPixelFormat().SetBitsAllocated(change_image.GetPixelFormat().GetBitsAllocated());
    im->GetPixelFormat().SetBitsStored(change_image.GetPixelFormat().GetBitsStored());
    im->GetPixelFormat().SetHighBit(change_image.GetPixelFormat().GetHighBit());
    im->GetPixelFormat().SetPixelRepresentation(
        change_image.GetPixelFormat().GetPixelRepresentation());
    im->SetSlope(change_image.GetSlope());
    im->SetIntercept(change_image.GetIntercept());

    // gdcm::Image im = change_image;
    ds = reader.GetFile().GetDataSet();
    im->SetDataElement(pixeldata);
    gdcm::UIDGenerator uid;
    gdcm::Attribute<0x0008, 0x18> ss;
    ss.SetValue(uid.Generate());
    ds.Replace(ss.GetAsDataElement());

    gdcm::Attribute<0x0020, 0x000e> ss2;
    ss2.SetValue(uid.Generate());
    ds.Replace(ss2.GetAsDataElement());


      snprintf(outputfilename, 1024 - 1, "%s/with_text/%08d.dcm", output.c_str(), i);
      gdcm::ImageWriter writer;
      writer.SetImage(*im);
      writer.SetFile(reader.GetFile());
      // snprintf(outputfilename, 1024 - 1, "%s/with_text/%08d.dcm", output.c_str(), i);
      writer.SetFileName(outputfilename);
      if (!writer.Write()) {
        return 1;
      }
    delete[] buffer;
  }

  return 0;
}

/* EOF */