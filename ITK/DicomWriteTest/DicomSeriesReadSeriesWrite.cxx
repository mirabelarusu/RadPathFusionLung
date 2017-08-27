#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkMetaDataObject.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include <vector>
#include "itksys/SystemTools.hxx"
#include "itkMinimumMaximumImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

typedef unsigned short  PixelType;
const unsigned int      Dimension = 3;
typedef itk::Image< PixelType, Dimension >      ImageType;
typedef itk::ImageSeriesReader<ImageType >      ReaderType;
typedef itk::GDCMImageIO                        ImageIOType;
typedef itk::GDCMSeriesFileNames                NamesGeneratorType;

void writeCharacteristics(itk::Image<PixelType, Dimension> * img)
{

  std::cout << "Size:   " << img->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "Spacing:" << img->GetSpacing() << std::endl;

  typedef itk::MinimumMaximumImageFilter< ImageType > CalcType;
  typename CalcType::Pointer dicomCalc = CalcType::New();
  dicomCalc->SetInput(img);
  dicomCalc->Update();
 
  unsigned short minVol = dicomCalc->GetMinimum();
  unsigned short maxVol = dicomCalc->GetMaximum();
  std::cout << "Vol Min:" << minVol   << std::endl;
  std::cout << "Vol Max:" << maxVol   << std::endl;
}

ReaderType::Pointer readDicom(std::string dicomDir, 
  itk::Image<PixelType, Dimension>::Pointer outImg,
  itk::GDCMImageIO::Pointer gdcmIO, 
  NamesGeneratorType::Pointer namesGenerator, 
  bool verbose = false )
{
  namesGenerator->SetInputDirectory( dicomDir );

  const ReaderType::FileNamesContainer & filenames =
                            namesGenerator->GetInputFileNames();

  std::size_t numberOfFileNames = filenames.size();
  if (verbose)
  {
    std::cout << numberOfFileNames << std::endl;
    for(unsigned int fni = 0; fni < numberOfFileNames; ++fni)
      {
      std::cout << "filename # " << fni << " = ";
      std::cout << filenames[fni] << std::endl;
      }
  }

  ReaderType::Pointer reader = ReaderType::New();

  reader->SetImageIO( gdcmIO );
  reader->SetFileNames( filenames );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cerr << "Exception thrown while writing the image" << std::endl;
    std::cerr << excp << std::endl;
    return NULL;
    }
  outImg = reader->GetOutput();

  writeCharacteristics(outImg);

  return reader;
}

void writeTags (itk::GDCMImageIO::Pointer gdcmIO)
{

  typedef itk::MetaDataDictionary   DictionaryType;
  DictionaryType & dictionary = gdcmIO->GetMetaDataDictionary();
  typedef itk::MetaDataObject< std::string > MetaDataStringType;


  DictionaryType::ConstIterator itr = dictionary.Begin();
  DictionaryType::ConstIterator end = dictionary.End();

  while( itr != end )
    {
    itk::MetaDataObjectBase::Pointer  entry = itr->second;
    MetaDataStringType::Pointer entryvalue =
      dynamic_cast<MetaDataStringType *>( entry.GetPointer() );
    if( entryvalue )
      {
      std::string tagkey   = itr->first;
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      std::cout << tagkey <<  " = " << tagvalue << std::endl;
      }
    ++itr;
    }
}

int main( int argc, char* argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Usage: " << argv[0] <<
      " DicomDirectory InputMhaTesting " <<
      "InputMhaToConvertToDicom OutputDicomDirectory NewSubjectName " << std::endl;
    return EXIT_FAILURE;
    }


  typename ImageType::Pointer dicomVol       = ImageType::New();
  ImageIOType::Pointer gdcmIO                = ImageIOType::New();
  NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
  ReaderType::Pointer reader = readDicom( argv[1], dicomVol, gdcmIO, namesGenerator);

  writeTags(gdcmIO);

  const ReaderType::FileNamesContainer & filenames =
                            namesGenerator->GetInputFileNames();
  unsigned int nbSlices = filenames.size();
	ReaderType::DictionaryRawPointer dictionary;
	ReaderType::DictionaryArrayType outputArray;
  for (unsigned int i = 0; i < nbSlices; i++)
  {
    dictionary = (*(reader->GetMetaDataDictionaryArray()))[i];

    std::string patientName = argv[5];
    std::string entryId( "0010|0010" );
    std::string value( patientName );
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, patientName );

    entryId = "0010|0020" ;
    std::string value2( patientName );
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, patientName );
   
    //Exam type
    std::string method = argv[6];
    entryId = "0008|0060" ;
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, method );
   
    std::string seriesDescription = argv[7];
    entryId =  "0008|103e";
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, seriesDescription );

    entryId =  "0070|0081";
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, seriesDescription );
    
    entryId = "0018|0015";
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, method );


    std::string manufacturersModelName = "https://github.com/mirabelarusu/RadPathFusionLung";
    entryId = "0008|1090" ;
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, manufacturersModelName );


    std::string manufacturer = argv[8];
    entryId = "0008|0070" ;
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, manufacturersModelName );

    std::string seriesNo = argv[9];
    entryId = "0020|0011" ;
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, seriesNo);

    //acquisition number
    entryId = "0020|0012" ;
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, seriesNo);



    if (atoi(argv[10])==0) //don't keep date
    {
      std::string newDate = argv[11];
      entryId = "0008|0020" ;
      itk::EncapsulateMetaData<std::string>( *dictionary, entryId, newDate );
    }
    

    //delete tag content 
    //Software
    entryId = "0018|1020";
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, "" );

    entryId = "0018|1030";
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, "" );

    entryId = "0018|1210";
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, "" );

    entryId = "0018|5100";
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, "" );

    entryId = "0018|1160";
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, "" );

    entryId = "0018|1140";
    itk::EncapsulateMetaData<std::string>( *dictionary, entryId, "" );

    outputArray.push_back(dictionary);
  }  

  ReaderType::Pointer reader2 = ReaderType::New();

  //testing whether the read dicom has the write size compared to the exported 
  //mha used for processing
  const char * inFilenameTest = argv[2];
  reader2->SetFileName( inFilenameTest  );
  try
    {
    reader2->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cerr << "Exception thrown while writing the image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  typename ImageType::Pointer mhaVol = reader2->GetOutput();

  writeCharacteristics(mhaVol);

  ReaderType::Pointer reader3 = ReaderType::New();
  //testing whether the read dicom has the write size compared to the exported 
  //mha used for processing
  const char * inFilename = argv[3];
  reader3->SetFileName( inFilename );
  try
    {
    reader3->Update();
    }
  catch (itk::ExceptionObject &excp)
    {
    std::cerr << "Exception thrown while writing the image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  typename ImageType::Pointer volCrop = reader3->GetOutput();


  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageType;
  typename ResampleImageType::Pointer resampler = ResampleImageType::New();

  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double>
    InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

  resampler->SetInput(volCrop);
  resampler->SetSize(mhaVol->GetLargestPossibleRegion().GetSize());
  resampler->SetTransform(itk::IdentityTransform<double,Dimension>::New());
  resampler->SetInterpolator(interpolator);
  resampler->SetDefaultPixelValue(0);
  resampler->SetOutputSpacing(mhaVol->GetSpacing());
  resampler->SetOutputDirection(mhaVol->GetDirection());
  resampler->SetOutputOrigin(mhaVol->GetOrigin());
  resampler->Update();
  typename ImageType::Pointer vol = resampler->GetOutput();

  writeCharacteristics(vol);

  const char * outputDirectory = argv[4];
  itksys::SystemTools::MakeDirectory( outputDirectory );
  const unsigned int        OutputDimension = 2;

  typedef itk::Image< PixelType, OutputDimension >    Image2DType;
  typedef itk::ImageSeriesWriter<
                             ImageType, Image2DType >  SeriesWriterType;
  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();

  seriesWriter->SetInput( vol );
  seriesWriter->SetImageIO( gdcmIO );
  namesGenerator->SetOutputDirectory( outputDirectory );

  seriesWriter->SetFileNames( namesGenerator->GetOutputFileNames() );
  seriesWriter->SetMetaDataDictionaryArray( & outputArray);

  std::cout << reader->GetMetaDataDictionaryArray() << std::endl;
  try
    {
    seriesWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing the series " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }


  typename ImageType::Pointer dicomVol2       = ImageType::New();
  ImageIOType::Pointer gdcmIO2                = ImageIOType::New();
  NamesGeneratorType::Pointer namesGenerator2 = NamesGeneratorType::New();
  ReaderType::Pointer readerOut = readDicom( outputDirectory, dicomVol2, gdcmIO2, 
    namesGenerator2 );

  writeTags(gdcmIO2);

  return EXIT_SUCCESS;
}
