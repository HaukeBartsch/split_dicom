
/*
 ./split_dicom -d data -o /tmp/
*/

#include <filesystem>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

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
#include <gdcmStringFilter.h>

// boost libraries
#include "boost/date_time.hpp"
#include "boost/filesystem.hpp"
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "optionparser.h"

struct Arg : public option::Arg {
  static option::ArgStatus Required(const option::Option &option, bool) { return option.arg == 0 ? option::ARG_ILLEGAL : option::ARG_OK; }
  static option::ArgStatus Empty(const option::Option &option, bool) { return (option.arg == 0 || option.arg[0] == 0) ? option::ARG_OK : option::ARG_IGNORE; }
};

enum optionIndex { UNKNOWN, INPUT, OUTPUT, HELP };
const option::Descriptor usage[] = {{UNKNOWN, 0, "", "", option::Arg::None,
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
      if (extension == ".tar" || extension == ".gz" || extension == ".zip" || extension == ".tgz" || extension == ".bz2")
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
  std::string output = ""; // directory path
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
  std::srand(random_seed); // we don't need that one... 

  std::vector<std::string> files = listFilesSTD(input);
  //
  // lets read in a random DICOM image (check that it is a DICOM image...)
  //
  bool foundDICOM = false;
  int pickImageIdx;
  int maxIter = 10;

  for (int i = 0; i < files.size(); i++) { // while (maxIter > 0) {
    pickImageIdx = i;
    // to check for DICOM lets use the normal gdcm::Read and CanRead
    gdcm::Reader r;
    r.SetFileName(files[pickImageIdx].c_str());
    if (r.CanRead() != true) {
      maxIter--;
      fprintf(stdout, "read failed [%d] for \"%s\"... try next\n", maxIter, files[pickImageIdx].c_str());
      continue; // try again
    }

    gdcm::ImageReader reader;
    reader.SetFileName(files[pickImageIdx].c_str());
    fprintf(stdout, "Try to read \"%s\"\n", files[pickImageIdx].c_str());
    if (!reader.Read()) {
      maxIter--;
      fprintf(stdout, "read failed %s... try next\n", files[pickImageIdx].c_str());
      continue; // try again
    }
    fprintf(stdout, "  DICOM file: \"%s\"\n", files[pickImageIdx].c_str());

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

    int lenL = pixelsize * (xmax / 2) * ymax;
    int lenR = len - lenL;
    //fprintf(stdout, "%d %d %d\n", len, lenL, lenR);

    char *bufferL = new char[lenL];
    char *bufferR = new char[lenR];
    {
      signed short *bvals = (signed short *)buffer;
      signed short *bL = (signed short *)bufferL;
      signed short *bR = (signed short *)bufferR;
      for (int yi = 0; yi < ymax; yi++) {
        for (int xi = 0; xi < xmax; xi++) {
          int idx = yi * xmax + xi;
          int idx2 = yi * xmax/2;
          if (xi > xmax/2)
            idx2 += xi - (xmax / 2);
          else
            idx2 += xi;

          if (xi <= xmax / 2 ) {
            bL[idx2] = bvals[idx];
          } else {
            bR[idx2] = bvals[idx];
          }
        }
      }
    }

    // what is max and min here? (don't change by the overlay)
    float current_image_min_value, current_image_max_value;
    {
      signed short *bvals = (signed short *)buffer;
      current_image_min_value = (float)bvals[0];
      current_image_max_value = (float)bvals[0];
      for (int p = 1; p < xmax * ymax; p++) {
        if (current_image_min_value > bvals[p])
          current_image_min_value = bvals[p];
        if (current_image_max_value < bvals[p])
          current_image_max_value = bvals[p];
      }
    }

    //
    // write the image pair
    //
    std::string dn = output;
    struct stat buf;
    if (!(stat(dn.c_str(), &buf) == 0)) {
      mkdir(dn.c_str(), 0777);
    }

    dn = output + "/left/";
    // struct stat buf;
    if (!(stat(dn.c_str(), &buf) == 0)) {
      mkdir(dn.c_str(), 0777);
    }
    dn = output + "/right/";
    if (!(stat(dn.c_str(), &buf) == 0)) {
      mkdir(dn.c_str(), 0777);
    }

    char outputfilename[1024];

//    gdcm::ImageWriter writer2;
//    writer2.SetImage(change.GetOutput());
//    writer2.SetFile(reader.GetFile());
//    snprintf(outputfilename, 1024 - 1, "%s/left/%08d.dcm", output.c_str(), i);
//    writer2.SetFileName(outputfilename);
//    if (!writer2.Write()) {
//      return 1;
//    }

    // now we can add the bitmap to the original data and write again
    // change_image.SetBuffer(buffer);
    // gdcm::DataElement pixeldata = change_image.GetDataElement();
    gdcm::DataElement pixeldataL(gdcm::Tag(0x7fe0, 0x0010));
    gdcm::DataElement pixeldataR(gdcm::Tag(0x7fe0, 0x0010));
    pixeldataL.SetByteValue(bufferL, lenL);
    pixeldataR.SetByteValue(bufferR, lenR);

    gdcm::SmartPointer<gdcm::Image> imL = new gdcm::Image;
    gdcm::SmartPointer<gdcm::Image> imR = new gdcm::Image;
    imL->SetNumberOfDimensions(2);
    imL->SetDimension(0, xmax/2);
    imL->SetDimension(1, ymax);
    fprintf(stdout, "dimensions are: %d %d\n", xmax / 2, ymax);
    imL->SetPhotometricInterpretation(change_image.GetPhotometricInterpretation());
    imL->GetPixelFormat().SetSamplesPerPixel(1);
    imL->GetPixelFormat().SetBitsAllocated(change_image.GetPixelFormat().GetBitsAllocated());
    imL->GetPixelFormat().SetBitsStored(change_image.GetPixelFormat().GetBitsStored());
    imL->GetPixelFormat().SetHighBit(change_image.GetPixelFormat().GetHighBit());
    imL->GetPixelFormat().SetPixelRepresentation(change_image.GetPixelFormat().GetPixelRepresentation());
    imL->SetSlope(change_image.GetSlope());
    imL->SetOrigin(change_image.GetOrigin());
    imL->SetIntercept(change_image.GetIntercept());
    //imL->SetWindowCenter(0);

    // we need to change the Origin as well so that we can overlay with referenced images
    imR->SetNumberOfDimensions(2);
    imR->SetDimension(0, xmax - (xmax/2));
    imR->SetDimension(1, ymax);
    fprintf(stdout, "dimensions are: %d %d\n", xmax - (xmax/2), ymax);
    imR->SetPhotometricInterpretation(change_image.GetPhotometricInterpretation());
    imR->GetPixelFormat().SetSamplesPerPixel(1);
    imR->GetPixelFormat().SetBitsAllocated(change_image.GetPixelFormat().GetBitsAllocated());
    imR->GetPixelFormat().SetBitsStored(change_image.GetPixelFormat().GetBitsStored());
    imR->GetPixelFormat().SetHighBit(change_image.GetPixelFormat().GetHighBit());
    imR->GetPixelFormat().SetPixelRepresentation(change_image.GetPixelFormat().GetPixelRepresentation());
    imR->SetSlope(change_image.GetSlope());
    imR->SetIntercept(change_image.GetIntercept());
    imR->SetOrigin(change_image.GetOrigin());

    // gdcm::Image im = change_image;
    ds = reader.GetFile().GetDataSet();
    imL->SetDataElement(pixeldataL);
    imR->SetDataElement(pixeldataR);
 
    // we should write everything new StudyInstanceUID, SeriesInstanceUID and SOPInstanceUID
 
    std::string StudyInstanceUID = "";
    gdcm::Tag StudyInstanceUID_tag(0x0020, 0x000d);
    if (ds.FindDataElement(StudyInstanceUID_tag)) {
      gdcm::StringFilter sf;
      sf.SetFile(file);
      StudyInstanceUID = sf.ToString(StudyInstanceUID_tag);
      //fprintf(stdout, "FOUND StudyInstanceUID: %s\n", StudyInstanceUID.c_str());
    }

    std::string SeriesInstanceUID = "";
    gdcm::Tag SeriesInstanceUID_tag(0x0020, 0x000e);
    if (ds.FindDataElement(SeriesInstanceUID_tag)) {
      gdcm::StringFilter sf;
      sf.SetFile(file);
      SeriesInstanceUID = sf.ToString(SeriesInstanceUID_tag);
      //fprintf(stdout, "FOUND StudyInstanceUID: %s\n", StudyInstanceUID.c_str());
    }

    std::string SeriesDescription = "";
    gdcm::Tag SeriesDescription_tag(0x0008, 0x103e);
    if (ds.FindDataElement(SeriesDescription_tag)) {
      gdcm::StringFilter sf;
      sf.SetFile(file);
      SeriesDescription = sf.ToString(SeriesDescription_tag);
      fprintf(stdout, "FOUND SeriesDescription: %s\n", SeriesDescription.c_str());
    }

    gdcm::Attribute<0x0020, 0x00d> ss3; // StudyInstanceUID - make this deterministic, use input + .1 and .2
    ss3.SetValue(StudyInstanceUID + ".1");
    ds.Replace(ss3.GetAsDataElement());

    gdcm::UIDGenerator uid;
    gdcm::Attribute<0x0008, 0x0018> ss;
    ss.SetValue(uid.Generate());
    gdcm::Tag SOPInstanceUID_tag(0x0008, 0x0018);
    if (ds.FindDataElement(SOPInstanceUID_tag)) {
      ds.Replace(ss.GetAsDataElement());
    }

    gdcm::Attribute<0x0020, 0x000e> ss2;
    ss2.SetValue(SeriesInstanceUID + ".1"); // same series instance uid for each file
    ds.Replace(ss2.GetAsDataElement());

    gdcm::Attribute<0x0008, 0x103e> ss6;
    ss6.SetValue(SeriesDescription + " (L)"); // same series instance uid for each file
    ds.Replace(ss6.GetAsDataElement());

    //gdcm::Attribute<0x28,0x8> at;
    //at.SetFromDataElement( ds.GetDataElement( at.GetTag() ) );
    //at.SetValue( at.GetValue() * 2 );
    //ds.Replace( at.GetAsDataElement() );

    snprintf(outputfilename, 1024 - 1, "%s/left/%08d.dcm", output.c_str(), i);
    gdcm::ImageWriter writer;
    writer.SetImage(*imL);
    writer.SetFile(reader.GetFile());
    writer.SetFileName(outputfilename);
    if (!writer.Write()) {
      return 1;
    }

    // now write right
    // gdcm::Attribute<0x0020, 0x00d> ss3; // StudyInstanceUID - make this deterministic, use input + .1 and .2
    ss3.SetValue(StudyInstanceUID + ".2");
    ds.Replace(ss3.GetAsDataElement());

    gdcm::Attribute<0x0008, 0x18> ss4;
    ss4.SetValue(uid.Generate());
    ds.Replace(ss4.GetAsDataElement());

    gdcm::Attribute<0x0020, 0x000e> ss5;
    ss5.SetValue(SeriesInstanceUID + ".2");
    ds.Replace(ss5.GetAsDataElement());

    gdcm::Attribute<0x0008, 0x103e> ss7;
    ss7.SetValue(SeriesDescription + " (R)"); // same series instance uid for each file
    ds.Replace(ss7.GetAsDataElement());

    snprintf(outputfilename, 1024 - 1, "%s/right/%08d.dcm", output.c_str(), i);
    gdcm::ImageWriter writer2;
    writer2.SetImage(*imR);
    writer2.SetFile(reader.GetFile());
    writer2.SetFileName(outputfilename);
    if (!writer2.Write()) {
      return 1;
    }


    delete[] buffer;
    delete[] bufferL;
    delete[] bufferR;
  }

  return 0;
}

/* EOF */