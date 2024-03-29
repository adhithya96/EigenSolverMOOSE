# MultiAppScalarToAuxScalarTransfer

!syntax description /Transfers/MultiAppScalarToAuxScalarTransfer

## Siblings transfer behavior

This transfer supports sending data from a MultiApp to a MultiApp if and only if the number of subapps
in the source MultiApp matches the number of subapps in the target MultiApp, and they are distributed
the same way on the parallel processes. Each source app is then matched to the target app with the same
subapp index.

## Example Input File Syntax

The following examples demonstrate the use the MultiAppScalarToAuxScalarTransfer for transferring data
to ([tosub]) and from ([fromsub]) sub-applications.

!listing multiapp_scalar_to_auxscalar_transfer/to_sub/parent.i block=Transfers id=tosub caption=Example use of MultiAppScalarToAuxScalarTransfer for transferring data +to+ sub-applications.

!listing multiapp_scalar_to_auxscalar_transfer/from_sub/parent.i block=Transfers id=fromsub caption=Example use of MultiAppScalarToAuxScalarTransfer for transferring data +from+ sub-applications.

!syntax parameters /Transfers/MultiAppScalarToAuxScalarTransfer

!syntax inputs /Transfers/MultiAppScalarToAuxScalarTransfer

!syntax children /Transfers/MultiAppScalarToAuxScalarTransfer
